// app/api/mle/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type MleRow = {
  root_prob_final: string | number[] | null;
  trans_prob_final: string | number[][] | number[] | null;
  root_name: string | null;
  d_pi_1: number | null;
  d_pi_2: number | null;
  d_pi_3: number | null;
  d_pi_4: number | null;
  d_m_1: number | null;
  d_m_2: number | null;
  d_m_3: number | null;
  d_m_4: number | null;
};

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || "").trim();
    const method = (searchParams.get("method") || "dirichlet").trim();

    if (!jobId) {
      return NextResponse.json({ error: "Missing job_id" }, { status: 400 });
    }

    const sql = `
      SELECT
        br.root_prob_final,
        br.trans_prob_final,
        br.root_name,
        j.d_pi_1, j.d_pi_2, j.d_pi_3, j.d_pi_4,
        j.d_m_1,  j.d_m_2,  j.d_m_3,  j.d_m_4
      FROM emtr_best_rep AS br
      LEFT JOIN emtr_jobs AS j ON j.job_id = br.job_id
      WHERE br.job_id = ? AND br.method = ?
      ORDER BY br.created_at DESC
      LIMIT 1
    `;

    const conn = db();
    const { rows } = await conn.execute<MleRow>(sql, [jobId, method]);

    if (!rows || rows.length === 0) {
      return NextResponse.json({
        job_id: jobId,
        method,
        root: null,
        trans: null,
        root_name: null,
        D_pi: null,
        D_M: null,
      });
    }

    const row = rows[0];

    // Parse root + trans which may be JSON strings
    const root =
      typeof row.root_prob_final === "string"
        ? (JSON.parse(row.root_prob_final) as number[])
        : row.root_prob_final;

    const transRaw =
      typeof row.trans_prob_final === "string"
        ? (JSON.parse(row.trans_prob_final) as unknown)
        : row.trans_prob_final;

    // Assemble Dirichlet alpha arrays
    const piCandidates = [row.d_pi_1, row.d_pi_2, row.d_pi_3, row.d_pi_4];
    const D_pi =
      piCandidates.every(v => typeof v === "number" && Number.isFinite(v) && v > 0)
        ? (piCandidates as number[])
        : null;

    const mCandidates = [row.d_m_1, row.d_m_2, row.d_m_3, row.d_m_4];
    const D_M =
      mCandidates.every(v => typeof v === "number" && Number.isFinite(v) && v > 0)
        ? (mCandidates as number[])
        : null;

    return NextResponse.json({
      job_id: jobId,
      method,
      root: Array.isArray(root) ? root : null,
      trans: transRaw ?? null,      // could be 4x4, flat-16, or map of 4x4s
      root_name: row.root_name ?? null,
      D_pi,                        
      D_M,                         
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ error: msg }, { status: 500 });
  }
}
