// app/api/jobs/[jobId]/rows/route.ts
import { NextResponse } from "next/server";
import mysql, { ResultSetHeader } from "mysql2/promise";

type Method = "main" | "dirichlet" | "parsimony" | "ssh";

export interface EmtrRow {
  method?: Method; // optional if you store method in jobs instead
  root: string;
  rep: number;
  iter: number;
  ll_pars: number;
  edc_ll_first: number;
  edc_ll_final: number;
  ll_final: number;
}

interface RowsPayload {
  rows: EmtrRow[];
}

/* ---------- tiny runtime validators (no any) ---------- */
function isFiniteNumber(x: unknown): x is number {
  return typeof x === "number" && Number.isFinite(x);
}

function isEmtrRow(obj: unknown): obj is EmtrRow {
  if (typeof obj !== "object" || obj === null) return false;
  const r = obj as Record<string, unknown>;
  const method = r.method;
  const methodOk =
    method === undefined ||
    method === "main" ||
    method === "dirichlet" ||
    method === "parsimony" ||
    method === "ssh";

  return (
    methodOk &&
    typeof r.root === "string" &&
    isFiniteNumber(r.rep) &&
    Number.isInteger(r.rep as number) &&
    isFiniteNumber(r.iter) &&
    Number.isInteger(r.iter as number) &&
    isFiniteNumber(r.ll_pars) &&
    isFiniteNumber(r.edc_ll_first) &&
    isFiniteNumber(r.edc_ll_final) &&
    isFiniteNumber(r.ll_final)
  );
}

function hasRowsProp(x: unknown): x is { rows: unknown } {
  return typeof x === "object" && x !== null && "rows" in (x as object);
}

function isRowsPayload(x: unknown): x is RowsPayload {
  if (!hasRowsProp(x)) return false;
  const rows = (x as { rows: unknown }).rows;
  return Array.isArray(rows) && rows.every(isEmtrRow);
}

/* ---------- DB connection (simple; you can switch to a pool) ---------- */
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

/* ---------- route handler ---------- */
export async function POST(
  req: Request,
  context: { params: Record<string, string | string[]> }
) {
  const raw = context.params["jobId"];
  const jobId = Array.isArray(raw) ? raw[0] : raw;
  if (typeof jobId !== "string" || jobId.length === 0) {
    return NextResponse.json(
      { ok: false, error: "Missing jobId" },
      { status: 400 }
    );
  }

  let bodyUnknown: unknown;
  try {
    bodyUnknown = await req.json();
  } catch {
    return NextResponse.json(
      { ok: false, error: "Invalid JSON body" },
      { status: 400 }
    );
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

    // If you also add a `method` column in emtr_rows, include it here.
    const sql = `
      INSERT INTO emtr_rows
      (job_id, root, rep, iter, ll_pars, edc_ll_first, edc_ll_final, ll_final)
      VALUES (?, ?, ?, ?, ?, ?, ?, ?)
      ON DUPLICATE KEY UPDATE job_id = job_id
    `;

    for (const r of rows) {
      await conn.execute<ResultSetHeader>(sql, [
        jobId,
        r.root,
        r.rep,
        r.iter,
        r.ll_pars,
        r.edc_ll_first,
        r.edc_ll_final,
        r.ll_final,
      ]);
    }

    await conn.commit();
    return NextResponse.json({ ok: true, inserted: rows.length });
  } catch (e: unknown) {
    await conn.rollback();
    const msg =
      typeof e === "object" && e && "message" in e
        ? String((e as { message?: unknown }).message)
        : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  } finally {
    await conn.end();
  }
}
