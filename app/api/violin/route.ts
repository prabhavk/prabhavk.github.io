// app/api/violin/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type MethodName = "Parsimony" | "Dirichlet" | "HSS";
type DBRow = { method: string; val: number };
type RootRow = { root: string | null };

function normalizeMethod(s: string): MethodName | null {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "Parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "Dirichlet";
  if (m.includes("hss") || m === "ssh") return "HSS";
  return null;
}

type ApiResp = {
  job_id: string;
  roots: string[];
  root?: string;
  series: Record<MethodName, number[]>;
  counts: Record<MethodName, number>;
  metric: "final" | "init";
};

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job = (searchParams.get("job") ?? "").trim();
    const rootParam = (searchParams.get("root") ?? "").trim();
    const metricParam = (searchParams.get("metric") ?? "final").trim().toLowerCase();
    const metric: "final" | "init" = metricParam === "init" ? "init" : "final";

    if (!job) {
      return NextResponse.json({ error: "Missing ?job=<job_id>" }, { status: 400 });
    }

    // 1) distinct roots for this job (completed runs, value present in the chosen metric)
    const col = metric === "init" ? "ll_init" : "ll_final";
    const rootsRows = await query<RootRow>(
      `
      SELECT DISTINCT r.root
        FROM emtr_init_final r
        JOIN emtr_jobs j ON j.job_id = r.job_id
       WHERE r.job_id = ?
         AND r.root IS NOT NULL
         AND r.root <> ''
         AND r.${col} IS NOT NULL
         AND j.status = 'completed'
       ORDER BY r.root
      `,
      [job]
    );
    const roots = rootsRows.map(r => r.root).filter((v): v is string => !!v && v.length > 0);

    // pick root: use client’s if valid, else first (if any)
    let root: string | undefined = undefined;
    if (rootParam && roots.includes(rootParam)) root = rootParam;
    else if (roots.length > 0) root = roots[0];

    // 2) series (filtered by chosen root if provided)
    const params: Array<string> = [job];
    const rootFilterSQL = root ? "AND r.root = ?" : "";
    if (root) params.push(root);

    const rows = await query<DBRow>(
      `
      SELECT r.method, CAST(r.${col} AS DOUBLE) AS val
        FROM emtr_init_final r
        JOIN emtr_jobs j ON j.job_id = r.job_id
       WHERE r.job_id = ?
         ${rootFilterSQL}
         AND r.${col} IS NOT NULL
         AND j.status = 'completed'
      `,
      params
    );

    const series: Record<MethodName, number[]> = { Parsimony: [], Dirichlet: [], HSS: [] };
    for (const r of rows) {
      const m = normalizeMethod(r.method);
      const v = Number(r.val);
      if (m && Number.isFinite(v)) series[m].push(v);
    }

    const counts = {
      Parsimony: series.Parsimony.length,
      Dirichlet: series.Dirichlet.length,
      HSS: series.HSS.length,
    };

    const payload: ApiResp = { job_id: job, roots, root, series, counts, metric };
    return NextResponse.json(payload);
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Internal error";
    return NextResponse.json({ error: message }, { status: 500 });
  }
}
