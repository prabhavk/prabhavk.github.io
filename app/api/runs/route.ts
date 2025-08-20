// app/api/runs/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query, type SQLParam } from "@/lib/db";

export const runtime = "nodejs";

const SORT_WHITELIST = new Set([
  "finished_at",
  "created_at",
  "dur_minutes",
  "reps",
  "thr",
  "max_iter",
  "job_id",
]);

type CountRow = { cnt: number };

// If you want, you can export this type and reuse on the client
export type RunRow = {
  job_id: string;
  status: string;
  created_at: string;
  finished_at: string;
  dur_minutes: number | null;
  thr: number | null;
  reps: number | null;
  max_iter: number | null;
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
    const page = Math.max(1, Number(searchParams.get("page") ?? 1));
    const pageSize = Math.min(100, Math.max(5, Number(searchParams.get("pageSize") ?? 20)));
    const sortParts = (searchParams.get("sort") ?? "finished_at:desc").split(":");
    const sortCol = SORT_WHITELIST.has(sortParts[0] ?? "") ? (sortParts[0] as string) : "finished_at";
    const sortDir = sortParts[1] === "asc" ? "ASC" : "DESC";

    const q = (searchParams.get("q") ?? "").trim();
    const from = searchParams.get("from"); // YYYY-MM-DD
    const to   = searchParams.get("to");   // YYYY-MM-DD

    const where: string[] = ["status = 'completed'"];
    const params: SQLParam[] = [];

    if (q) {
      where.push("job_id LIKE ?");
      params.push(`%${q}%`);
    }
    if (from) {
      where.push("finished_at >= ?");
      params.push(from);
    }
    if (to) {
      where.push("finished_at < DATE_ADD(?, INTERVAL 1 DAY)");
      params.push(to);
    }

    const whereSQL = where.length ? `WHERE ${where.join(" AND ")}` : "";

    // count
    const [{ cnt }] = await query<CountRow>(`SELECT COUNT(*) AS cnt FROM emtr_jobs ${whereSQL}`, params);
    const total = Number(cnt ?? 0);

    // page query (columns aligned to your schema sample)
    const pageParams: SQLParam[] = [...params, pageSize, (page - 1) * pageSize];

    const rows = await query<RunRow>(
      `SELECT
         job_id, status, created_at, finished_at, dur_minutes, thr, reps, max_iter,
         d_pi_1, d_pi_2, d_pi_3, d_pi_4,
         d_m_1, d_m_2, d_m_3, d_m_4
       FROM emtr_jobs
       ${whereSQL}
       ORDER BY ${sortCol} ${sortDir}
       LIMIT ? OFFSET ?`,
      pageParams
    );

    return NextResponse.json({ page, pageSize, total, rows });
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Internal error";
    console.error("/api/runs error", err);
    return NextResponse.json({ error: message }, { status: 500 });
  }
}
