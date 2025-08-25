// app/api/ecdll/reps/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Row = { rep: number | null };

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || "").trim();
    // optional; default to 'dirichlet' if you want to scope reps per method
    const method = (searchParams.get("method") || "").trim().toLowerCase();

    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }

    const conn = db();

    // Build SQL dynamically only to include method filter if provided
    let sql = `
      SELECT DISTINCT rep
      FROM emtr_all_info
      WHERE job_id = ?
        AND rep IS NOT NULL
    `;
    const params: Array<string | number> = [jobId];

    if (method) {
      sql += ` AND method = ?`;
      params.push(method);
    }

    sql += ` ORDER BY rep ASC`;

    const { rows } = await conn.execute<Row>(sql, params);

    const reps = (rows ?? [])
      .map((r) => (typeof r.rep === "number" && Number.isFinite(r.rep) ? r.rep : null))
      .filter((v): v is number => v !== null);

    return NextResponse.json({
      ok: true,
      job_id: jobId,
      method: method || null,
      reps,
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
