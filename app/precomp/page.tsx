// app/precomp/page.tsx
import { query as q } from "@/lib/db";
import PrecompClient from "./precompclient";

export const dynamic = "force-dynamic";

type JobRow = {
  job_id: string;
  created_at: string;
  status: string | null;
};

export default async function PrecompPage() {
  // Try canonical jobs table first; fall back to distinct job_ids in trees
  let jobs: JobRow[] = [];
  try {
    jobs = await q<JobRow>(
      `SELECT job_id, created_at, status
         FROM emtr_jobs
         ORDER BY created_at DESC
         LIMIT 200`
    );
  } catch {
    jobs = await q<JobRow>(
      `SELECT job_id,
              MAX(created_at) AS created_at,
              NULL AS status
         FROM emtr_trees
         GROUP BY job_id
         ORDER BY created_at DESC
         LIMIT 200`
    );
  }

  return <PrecompClient jobs={jobs} />;
}
