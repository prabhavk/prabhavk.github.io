// app/api/jobs/[jobId]/rows/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

// Row shape FROM THE WORKER
interface EmtrRowIn {
  method: string;                  // e.g., "Parsimony" | "Dirichlet" | "SSH" | "main"
  root: string;
  rep: number | string;
  iter: number | string;
  ll_init: number | string;
  ecd_ll_first: number | string;   // NOTE: corrected "ecd_*" names
  ecd_ll_final: number | string;
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

// Normalize worker-provided method to your DB enum (lowercase)
function normMethod(m: unknown): "parsimony" | "dirichlet" | "ssh" | "main" {
  const s = String(m ?? "").trim().toLowerCase();
  if (s === "parsimony" || s === "dirichlet" || s === "ssh" || s === "main") return s;
  // Accept common capitalizations from native logs:
  if (s === "reparameterized mle") return "ssh"; // if you ever emit that label
  return "main";
}

function isRow(x: unknown): x is EmtrRowIn {
  if (typeof x !== "object" || x === null) return false;
  const r = x as Record<string, unknown>;
  // Validate types and numeric fields
  return (
    typeof r.method !== "undefined" &&
    typeof r.root === "string" &&
    isSafeInt(r.rep) &&
    isSafeInt(r.iter) &&
    isFiniteNum(r.ll_init) &&
    isFiniteNum(r.ecd_ll_first) &&
    isFiniteNum(r.ecd_ll_final) &&
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
    const ct = req.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      return NextResponse.json({ ok: false, error: "Content-Type must be application/json" }, { status: 415 });
    }

    const { jobId } = await params;
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }

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

    // Validate + collect
    const valid: EmtrRowIn[] = [];
    for (const r of rows) {
      if (!isRow(r)) {
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
        r.method,       
        r.root,
        toInt(r.rep),
        toInt(r.iter),
        toNum(r.ll_init),
        toNum(r.ecd_ll_first),
        toNum(r.ecd_ll_final),
        toNum(r.ll_final),
      ]);

      // Columns reflect corrected "ecd_*" names
      const sql = `
        INSERT INTO emtr_rows
          (job_id, method, root, rep, iter, ll_init, ecd_ll_first, ecd_ll_final, ll_final)
        VALUES ${placeholders}
        ON DUPLICATE KEY UPDATE job_id = job_id
      `;

      await conn.execute(sql, args);
      inserted += batch.length;
    }

    return NextResponse.json({ ok: true, inserted, skipped: rows.length - inserted, total: rows.length });
  } catch (e) {
    const msg =
      e && typeof e === "object" && "message" in e ? String((e as { message?: unknown }).message) : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
