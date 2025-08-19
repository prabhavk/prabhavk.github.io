// app/api/jobs/[jobId]/rows/route.ts
import { NextResponse } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Method = "main" | "dirichlet" | "parsimony" | "ssh";

interface EmtrRow {
  method: Method;
  root: string;
  rep: number;        // integer >= 0
  iter: number;       // integer >= 0
  ll_pars: number;
  edc_ll_first: number;
  edc_ll_final: number;
  ll_final: number;
}

/* ---------- config ---------- */
const MAX_ROWS_PER_REQUEST = 5_000;
const INSERT_CHUNK_SIZE = 500;

/* ---------- utils ---------- */
const isSafeInt = (x: unknown) => typeof x === "number" && Number.isInteger(x) && x >= 0;
const isFiniteNum = (x: unknown) => typeof x === "number" && Number.isFinite(x);

function isEmtrRow(x: unknown): x is EmtrRow {
  if (typeof x !== "object" || x === null) return false;
  const r = x as Record<string, unknown>;
  const m = r.method;
  const methodOk = m === "main" || m === "dirichlet" || m === "parsimony" || m === "ssh";
  return (
    methodOk &&
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
  req: Request,
  { params }: { params: { jobId: string } }   // ✅ accept params from Next.js
) {
  try {
    // Require JSON
    const ct = req.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      return NextResponse.json({ ok: false, error: "Content-Type must be application/json" }, { status: 415 });
    }

    // jobId from route params
    const jobId = params.jobId;
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }

    // Parse body
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
    const valid: EmtrRow[] = [];
    for (const r of rows) {
      if (isEmtrRow(r)) valid.push(r);
      else return NextResponse.json({ ok: false, error: "Row validation failed" }, { status: 400 });
    }

    // Build SQL and insert in chunks
    const conn = db();
    let inserted = 0;

    for (const batch of chunk(valid, INSERT_CHUNK_SIZE)) {
      const valuesSql = batch.map(() => "(?,?,?,?,?,?,?,?,?)").join(",");
      const sqlParams = batch.flatMap((r) => [   // ✅ renamed from `params` → `sqlParams`
        jobId,
        r.method,
        r.root,
        r.rep,
        r.iter,
        r.ll_pars,
        r.edc_ll_first,
        r.edc_ll_final,
        r.ll_final,
      ]);

      const sql = `
        INSERT INTO emtr_rows
          (job_id, method, root, rep, iter, ll_pars, edc_ll_first, edc_ll_final, ll_final)
        VALUES ${valuesSql}
        ON DUPLICATE KEY UPDATE job_id = job_id
      `;

      await conn.execute(sql, sqlParams);
      inserted += batch.length;
    }

    return NextResponse.json({ ok: true, inserted, skipped: rows.length - inserted, total: rows.length });
  } catch (e) {
    const msg =
      e && typeof e === "object" && "message" in e ? String((e as { message?: unknown }).message) : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
