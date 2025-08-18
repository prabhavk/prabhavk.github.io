import { NextRequest, NextResponse } from "next/server";
import { pool } from "@/lib/db";

export const runtime = "nodejs";
export const dynamic = "force-dynamic";

type Row = {
  init_method: "Parsimony" | "Dirichlet" | "SSH";
  root: string;
  repetition: number;
  iter: number;
  ll_initial: number;
  ecd_ll_first: number;
  ecd_ll_final: number;
  ll_final: number;
};

export async function POST(
  req: NextRequest,
  ctx: { params: Promise<{ jobId: string }> }  
) {
  const { jobId } = await ctx.params;          

  let payload: { rows: Row[] };
  try {
    payload = await req.json();
  } catch {
    return NextResponse.json({ ok: false, error: "invalid json" }, { status: 400 });
  }

  const rows = payload?.rows;
  if (!Array.isArray(rows) || rows.length === 0) {
    return NextResponse.json({ ok: false, error: "rows[] required" }, { status: 400 });
  }

  const conn = await pool.getConnection();
  try {
    await conn.beginTransaction();

    const sql = `
      INSERT INTO emtr_rows
        (job_id, init_method, root, repetition, iter, ll_initial, ecd_ll_first, ecd_ll_final, ll_final)
      VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
      ON DUPLICATE KEY UPDATE job_id = job_id
    `;

    for (const r of rows) {
      await conn.execute(sql, [
        jobId, r.init_method, r.root, r.repetition, r.iter,
        r.ll_initial, r.ecd_ll_first, r.ecd_ll_final, r.ll_final,
      ]);
    }

    await conn.commit();
    return NextResponse.json({ ok: true, inserted: rows.length });
  } catch (e: any) {
    await conn.rollback();
    return NextResponse.json({ ok: false, error: String(e) }, { status: 500 });
  } finally {
    conn.release();
  }
}
