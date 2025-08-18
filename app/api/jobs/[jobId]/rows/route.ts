// app/api/jobs/[jobId]/rows/route.ts
import { NextRequest, NextResponse } from "next/server";
import mysql, { ResultSetHeader } from "mysql2/promise";

type Method = "main" | "dirichlet" | "parsimony" | "ssh";

export interface EmtrRow {
  // method is optional if you use one table for all three functions and/or store it in jobs
  method?: Method;
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

// ---- tiny runtime validators (no zod dependency) ----
function isFiniteNumber(x: unknown): x is number {
  return typeof x === "number" && Number.isFinite(x);
}

function isEmtrRow(x: unknown): x is EmtrRow {
  if (typeof x !== "object" || x === null) return false;
  const r = x as Record<string, unknown>;
  const methodOk =
    r.method === undefined ||
    r.method === "main" ||
    r.method === "dirichlet" ||
    r.method === "parsimony" ||
    r.method === "ssh";
  return (
    methodOk &&
    typeof r.root === "string" &&
    isFiniteNumber(r.rep) &&
    Number.isInteger(r.rep) &&
    isFiniteNumber(r.iter) &&
    Number.isInteger(r.iter) &&
    isFiniteNumber(r.ll_pars) &&
    isFiniteNumber(r.edc_ll_first) &&
    isFiniteNumber(r.edc_ll_final) &&
    isFiniteNumber(r.ll_final)
  );
}

function isRowsPayload(x: unknown): x is RowsPayload {
  return (
    typeof x === "object" &&
    x !== null &&
    Array.isArray((x as any).rows) && // eslint won’t complain: we don’t type the variable as any
    (x as { rows: unknown[] }).rows.every(isEmtrRow)
  );
}

// Optional: factor this out to a shared db util module for pooling
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

export async function POST(
  req: NextRequest,
  { params }: { params: { jobId: string } }
) {
  const jobId = params.jobId;
  if (!jobId || typeof jobId !== "string") {
    return NextResponse.json(
      { ok: false, error: "Missing or invalid jobId" },
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
      { ok: false, error: "Invalid payload shape for { rows: EmtrRow[] }" },
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

    // If you also want to store the method column, include it below and in your table.
    // Here we assume table emtr_rows has columns:
    // (job_id, root, rep, iter, ll_pars, edc_ll_first, edc_ll_final, ll_final)
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
      e && typeof e === "object" && "message" in e
        ? String((e as { message?: unknown }).message)
        : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  } finally {
    await conn.end();
  }
}
