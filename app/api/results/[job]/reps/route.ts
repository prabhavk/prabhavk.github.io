// app/api/results/[job]/reps/route.ts
import { NextResponse } from "next/server";
import { query as q } from "@/lib/db";

export const runtime = "nodejs";

type Row = { rep: number | string };

export async function GET(req: Request) {
  const { pathname } = new URL(req.url);
  // Expect: /api/results/<job>/reps
  const m = pathname.match(/\/api\/results\/([^/]+)\/reps\/?$/);
  const job = decodeURIComponent(m?.[1] ?? "");

  if (!job) {
    return NextResponse.json({ ok: false, job_id: null, reps: [] });
  }

  const rows = await q<Row>(
    `SELECT DISTINCT rep
       FROM emtr_trifle_layers
      WHERE job_id = ?
      ORDER BY rep ASC`,
    [job]
  );

  const reps = (rows ?? []).map(r => Number(r.rep)).filter(Number.isFinite);
  return NextResponse.json({ ok: true, job_id: job, reps });
}
