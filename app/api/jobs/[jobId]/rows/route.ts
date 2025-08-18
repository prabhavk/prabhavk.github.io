// app/api/jobs/[jobId]/rows/route.ts
import { NextResponse } from "next/server";
import mysql, { ResultSetHeader } from "mysql2/promise";

type Method = "main" | "dirichlet" | "parsimony" | "ssh";

export interface EmtrRow {
  method?: Method;
  root: string;
  rep: number;
  iter: number;
  ll_pars: number;
  edc_ll_first: number;
  edc_ll_final: number;
  ll_final: number;
}

/* ---------- tiny runtime validators (no any) ---------- */
function isFiniteNumber(x: unknown): x is number {
  return typeof x === "number" && Number.isFinite(x);
}
function isEmtrRow(obj: unknown): obj is EmtrRow {
  if (typeof obj !== "object" || obj === null) return false;
  const r = obj as Record<string, unknown>;
  const m = r.method;
  const methodOk =
    m === undefined || m === "main" || m === "dirichlet" || m === "parsimony" || m === "ssh";
  return (
    methodOk &&
    typeof r.root === "string" &&
    isFiniteNumber(r.rep) && Number.isInteger(r.rep as number) &&
    isFiniteNumber(r.iter) && Number.isInteger(r.iter as number) &&
    isFiniteNumber(r.ll_pars) &&
    isFiniteNumber(r.edc_ll_first) &&
    isFiniteNumber(r.edc_ll_final) &&
    isFiniteNumber(r.ll_final)
  );
}
function isRowsPayload(x: unknown): x is { rows: EmtrRow[] } {
  if (typeof x !== "object" || x === null) return false;
  const rows = (x as { rows?: unknown }).rows;
  return Array.isArray(rows) && rows.every(isEmtrRow);
}

/* ---------- DB conn ---------- */
async function getConn() {
  return mysql.createConnection({
    host: process.env.EMTR_DB_HOST ?? "127.0.0.1",
    user: process.env.EMTR_DB_USER ?? "emtr",
    password: process.env.EMTR_DB_PASS ?? "emtrpw",
    database: process.env.EMTR_DB_NAME ?? "emtr",
    port: Number(process.env.EMTR_DB_PORT ?? 3306),
    multipleStatements: false,
  });
}

/* ---------- handler (single-arg; parse jobId from URL) ---------- */
export async function POST(req: Request) {
  const { pathname } = new URL(req.url);
  // matches /api/jobs/{jobId}/rows  (with or without trailing slash)
  const m = pathname.match(/\/api\/jobs\/([^/]+)\/rows\/?$/);
  const jobId = m?.[1] ?? "";
  if (!jobId) {
    return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
  }

  let bodyUnknown: unknown;
  try {
    bodyUnknown = await req.json();
  } catch {
    return NextResponse.json({ ok: false, error: "Invalid JSON body" }, { status: 400 });
  }

  if (!isRowsPayload(bodyUnknown)) {
    return NextResponse.json(
      { ok: false, error: "Invalid payload; expected { rows: EmtrRow[] }" },
      { status: 400 }
    );
  }

  const { rows } = bodyUnknown;
  if (rows.length === 0) {
    return NextResponse.json({ ok: true, inserted: 0 });
  }

  const conn = await getConn();
  try {
    await conn.beginTransaction();
    const sql = `
      INSERT INTO emtr_rows
      (job_id, root, rep, iter, ll_pars, edc_ll_first, edc_ll_final, ll_final)
      VALUES (?, ?, ?, ?, ?, ?, ?, ?)
      ON DUPLICATE KEY UPDATE job_id = job_id
    `;
    for (const r of rows) {
      await conn.execute<ResultSetHeader>(sql, [
        jobId, r.root, r.rep, r.iter, r.ll_pars, r.edc_ll_first, r.edc_ll_final, r.ll_final,
      ]);
    }
    await conn.commit();
    return NextResponse.json({ ok: true, inserted: rows.length });
  } catch (e: unknown) {
    await conn.rollback();
    const msg = typeof e === "object" && e && "message" in e
      ? String((e as { message?: unknown }).message)
      : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  } finally {
    await conn.end();
  }
}
