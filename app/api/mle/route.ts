// app/api/mle/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type MleRow = {
  root_prob_final: string | number[] | null;
  trans_prob_final: string | number[][] | number[] | null;
  // include if you eventually select them:
  // root_name?: string | null;
  // job_id?: string;
  // method?: string;
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
      SELECT root_prob_final, trans_prob_final
      FROM emtr_best_rep
      WHERE job_id = ? AND method = ?
      ORDER BY created_at DESC
      LIMIT 1
    `;

    const conn = db();
    // PlanetScale client's execute<T> returns { rows: T[] }
    const { rows } = await conn.execute<MleRow>(sql, [jobId, method]);

    if (!rows || rows.length === 0) {
      return NextResponse.json({ root: null, trans: null, job_id: jobId, method });
    }

    const row = rows[0];

    const root =
      typeof row.root_prob_final === "string"
        ? (JSON.parse(row.root_prob_final) as number[])
        : row.root_prob_final;

    const transRaw =
      typeof row.trans_prob_final === "string"
        ? (JSON.parse(row.trans_prob_final) as unknown)
        : row.trans_prob_final;

    return NextResponse.json({
      job_id: jobId,
      method,
      root: Array.isArray(root) ? root : null,
      trans: transRaw ?? null, // can be 4x4 array, flat 16, or dict of 4x4s (handled in UI)
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ error: msg }, { status: 500 });
  }
}
