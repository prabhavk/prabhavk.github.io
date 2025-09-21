// app/api/results/[job]/reps/roots/route.ts
import { NextResponse } from "next/server";
import { query as q } from "@/lib/db";

export const runtime = "nodejs";

type Row = { root: string | null };

export async function GET(req: Request) {
  const { pathname, search } = new URL(req.url);
  // Expect path like: /api/results/<job>/reps/roots
  const m = pathname.match(/\/api\/results\/([^/]+)\/reps\/roots\/?$/);
  const job = decodeURIComponent(m?.[1] ?? "");

  const u = new URL(req.url);
  const rep = Number(u.searchParams.get("rep"));

  if (!job || !Number.isFinite(rep)) {
    return NextResponse.json({ ok: false, job_id: job || null, rep: Number.isFinite(rep) ? rep : null, roots: [] });
  }

  const rows = await q<Row>(
    `SELECT DISTINCT root
       FROM emtr_trifle_layers
      WHERE job_id = ? AND rep = ?
      ORDER BY root ASC`,
    [job, rep]
  );

  const roots = (rows ?? []).map(r => r.root).filter((s): s is string => !!s);
  return NextResponse.json({ ok: true, job_id: job, rep, roots });
}
