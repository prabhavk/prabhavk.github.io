// app/api/violin/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type MethodName = "Parsimony" | "Dirichlet" | "SSH";
type DBRow = { method: string; ll_final: number };
type RootRow = { root: string | null };

function normalizeMethod(s: string): MethodName | null {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "Parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "Dirichlet";
  if (m.includes("ssh")) return "SSH";
  return null;
}

type ApiResp = {
  job_id: string;
  roots: string[];                     // <-- new
  root?: string;                       // (server-selected if none provided)
  series: Record<MethodName, number[]>;
  counts: Record<MethodName, number>;
};

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job = (searchParams.get("job") ?? "").trim();
    const rootParam = (searchParams.get("root") ?? "").trim();

    if (!job) {
      return NextResponse.json({ error: "Missing ?job=<job_id>" }, { status: 400 });
    }

    // 1) Fetch distinct roots for this job (completed runs only)
    const rootsRows = await query<RootRow>(
      `
      SELECT DISTINCT r.root
        FROM emtr_init_final r
        JOIN emtr_jobs j ON j.job_id = r.job_id
       WHERE r.job_id = ?
         AND r.root IS NOT NULL
         AND j.status = 'completed'
       ORDER BY r.root
      `,
      [job]
    );
    const roots = rootsRows
      .map((r) => r.root)
      .filter((v): v is string => !!v && v.length > 0);

    // 2) Choose root:
    //    - if client provided one and it's valid, use it
    //    - else if any exist, default to the first
    //    - else leave undefined (means "no roots for this job")
    let root: string | undefined = undefined;
    if (rootParam && roots.includes(rootParam)) {
      root = rootParam;
    } else if (roots.length > 0) {
      root = roots[0];
    }

    // 3) Fetch series (filtered by chosen root if available)
    const params: Array<string> = [job];
    const rootFilterSQL = root ? "AND r.root = ?" : "";
    if (root) params.push(root);

    const rows = await query<DBRow>(
      `
      SELECT r.method, r.ll_final
        FROM emtr_init_final r
        JOIN emtr_jobs j ON j.job_id = r.job_id
       WHERE r.job_id = ?
         ${rootFilterSQL}
         AND r.ll_final IS NOT NULL
         AND j.status = 'completed'
      `,
      params
    );

    const series: ApiResp["series"] = { Parsimony: [], Dirichlet: [], SSH: [] };
    for (const r of rows) {
      const m = normalizeMethod(r.method);
      const v = Number(r.ll_final);
      if (m && Number.isFinite(v)) series[m].push(v);
    }

    const counts: ApiResp["counts"] = {
      Parsimony: series.Parsimony.length,
      Dirichlet: series.Dirichlet.length,
      SSH: series.SSH.length,
    };

    const payload: ApiResp = { job_id: job, roots, root, series, counts };
    return NextResponse.json(payload);
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Internal error";
    return NextResponse.json({ error: message }, { status: 500 });
  }
}
