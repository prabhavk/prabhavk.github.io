// app/api/results/route.ts
import { NextResponse } from "next/server";
import { query } from "@/lib/db";

type JobRow = { job_id: string };

export async function GET(_req: Request) {
  // Ask for the row *shape*; query returns JobRow[]
  const rows = await query<JobRow>(
    `
      SELECT DISTINCT job_id
        FROM emtr_trifle_runs
      UNION
      SELECT DISTINCT job_id
        FROM emtr_trifle_layers
      ORDER BY job_id DESC
    `
  );

  const jobs: string[] = rows
    .map((r) => r.job_id)
    .filter((id): id is string => typeof id === "string" && id.length > 0);

  return NextResponse.json({ jobs });
}
