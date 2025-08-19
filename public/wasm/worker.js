import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

// Row shape AFTER worker normalization (numbers may arrive as strings)
interface EmtrRowIn {
  root: string;
  rep: number | string;
  iter: number | string;
  ll_pars: number | string;
  edc_ll_first: number | string;
  edc_ll_final: number | string;
  ll_final: number | string;
}

/* ---------- config ---------- */
const MAX_ROWS_PER_REQUEST = 5_000;
const INSERT_CHUNK_SIZE = 500;

/* ---------- utils ---------- */
const toInt = (x: unknown) => (typeof x === "number" ? x : parseInt(String(x), 10));
const toNum = (x: unknown) => (typeof x === "number" ? x : Number(String(x)));
const isFiniteNum = (x: unknown) => Number.isFinite(toNum(x));
const isSafeInt = (x: unknown) => Number.isInteger(toInt(x)) && toInt(x) >= 0;

function isRow(x: unknown): x is EmtrRowIn {
  if (typeof x !== "object" || x === null) return false;
  const r = x as Record<string, unknown>;
  return (
    typeof r.root === "string" &&
    isSafeInt(r.rep) &&
    isSafeInt(r.iter) &&
    isFiniteNum(r.ll_pars) &&
    isFiniteNum(r.edc_ll_first) &&
    isFiniteNum(r.edc_ll_final) &&
    isFiniteNum(r.ll_final)
  );
}

function chunk<T>(arr: T[], size: number): T[][] {
  const out: T[][] = [];
  for (let i = 0; i < arr.length; i += size) out.push(arr.slice(i, i + size));
  return out;
}

/* ---------- handler ---------- */
export async function POST(
  req: NextRequest,
  // Next.js 15: params is a Promise
  { params }: { params: Promise<{ jobId: string }> }
) {
  try {
    // Content-Type guard
    const ct = req.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      return NextResponse.json({ ok: false, error: "Content-Type must be application/json" }, { status: 415 });
    }

    const { jobId } = await params;
    if (!jobId) return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });

    const body = await req.json().catch(() => null as unknown);
    const rows = (body && typeof body === "object" && Array.isArray((body as { rows?: unknown[] }).rows))
      ? (body as { rows: unknown[] }).rows
      : null;

    if (!rows) {
      return NextResponse.json({ ok: false, error: "Invalid payload; expected { rows: EmtrRow[] }" }, { status: 400 });
    }
    if (rows.length === 0) {
      return NextResponse.json({ ok: true, inserted: 0, skipped: 0, total: 0 });
    }
    if (rows.length > MAX_ROWS_PER_REQUEST) {
      return NextResponse.json(
        { ok: false, error: `Too many rows; max ${MAX_ROWS_PER_REQUEST} per request` },
        { status: 413 }
      );
    }

    // Validate
    const valid: EmtrRowIn[] = [];
    for (const r of rows) {
      if (!isRow(r)) {
        // Optional: include small sample of bad row for debugging
        return NextResponse.json({ ok: false, error: "Row validation failed" }, { status: 400 });
      }
      valid.push(r);
    }

    const conn = db();
    let inserted = 0;

    for (const batch of chunk(valid, INSERT_CHUNK_SIZE)) {
      const placeholders = batch.map(() => "(?,?,?,?,?,?,?,?,?)").join(",");
      const args = batch.flatMap((r) => [
        jobId,
        "main",                // force method to "main"
        r.root,
        toInt(r.rep),
        toInt(r.iter),
        toNum(r.ll_pars),
        toNum(r.edc_ll_first),
        toNum(r.edc_ll_final),
        toNum(r.ll_final),
      ]);

      // Use ON DUPLICATE KEY only if you have a UNIQUE key on (job_id, root, rep, iter)
      const sql = `
        INSERT INTO emtr_rows
          (job_id, method, root, rep, iter, ll_pars, edc_ll_first, edc_ll_final, ll_final)
        VALUES ${placeholders}
        ON DUPLICATE KEY UPDATE job_id = job_id
      `;
      await conn.execute(sql, args);
      inserted += batch.length;
    }

    return NextResponse.json({ ok: true, inserted, skipped: rows.length - inserted, total: rows.length });
  } catch (e: any) {
    return NextResponse.json({ ok: false, error: String(e?.message ?? e ?? "Unknown error") }, { status: 500 });
  }
}
