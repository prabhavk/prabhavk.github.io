// app/api/ecdll/reps/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Row = { rep: number | null };

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || "").trim();
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }

    const sql = `
      SELECT DISTINCT rep
      FROM emtr_best_rep
      WHERE job_id = ?
      ORDER BY rep ASC
    `;
    const conn = db();
    const { rows } = await conn.execute<Row>(sql, [jobId]);

    const reps = rows
      .map(r => (r.rep ?? undefined))
      .filter((v): v is number => typeof v === "number" && Number.isFinite(v));

    return NextResponse.json({ ok: true, job_id: jobId, reps });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
