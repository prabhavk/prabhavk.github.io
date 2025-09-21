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
  layer?: number;
};

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job = (searchParams.get("job") ?? "").trim();
    const rootParam = (searchParams.get("root") ?? "").trim();
    const metricParam = (searchParams.get("metric") ?? "final").trim().toLowerCase();
    const metric: "final" | "init" = metricParam === "init" ? "init" : "final";

    // NEW: layer (0,1,2), default 1 (medium)
    const layerParam = searchParams.get("layer");
    const layer = layerParam == null ? 1 : Math.max(0, Math.min(2, Number(layerParam) | 0));

    if (!job) {
      return NextResponse.json({ error: "Missing ?job=<job_id>" }, { status: 400 });
    }

    // Pick column names from emtr_trifle_layers
    const col = metric === "init" ? "ll_initial" : "ll_final";

    // 1) distinct roots for this job *and layer* with a non-null chosen metric
    //    (keeps behavior similar to your previous route, but now scoped by layer)
    const rootsRows = await query<RootRow>(
      `
      SELECT DISTINCT t.root
        FROM emtr_trifle_layers t
        JOIN emtr_jobs j ON j.job_id = t.job_id
       WHERE t.job_id = ?
         AND t.layer = ?
         AND t.root IS NOT NULL
         AND t.root <> ''
         AND t.${col} IS NOT NULL
         AND j.status = 'completed'
       ORDER BY t.root
      `,
      [job, layer]
    );
    const roots = (rootsRows ?? [])
      .map(r => r.root)
      .filter((v): v is string => !!v && v.length > 0);

    // pick root: use client’s if valid, else first (if any)
    let root: string | undefined = undefined;
    if (rootParam && roots.includes(rootParam)) root = rootParam;
    else if (roots.length > 0) root = roots[0];

    // 2) series (filtered by chosen root if provided), from emtr_trifle_layers
    const params: Array<string | number> = [job, layer];
    const rootFilterSQL = root ? "AND t.root = ?" : "";
    if (root) params.push(root);

    const rows = await query<DBRow>(
      `
      SELECT t.method, CAST(t.${col} AS DOUBLE) AS val
        FROM emtr_trifle_layers t
        JOIN emtr_jobs j ON j.job_id = t.job_id
       WHERE t.job_id = ?
         AND t.layer = ?
         ${rootFilterSQL}
         AND t.${col} IS NOT NULL
         AND j.status = 'completed'
       ORDER BY t.rep
      `,
      params
    );

    const series: Record<MethodName, number[]> = { Parsimony: [], Dirichlet: [], HSS: [] };
    for (const r of rows ?? []) {
      const m = normalizeMethod(r.method);
      const v = Number(r.val);
      if (m && Number.isFinite(v)) series[m].push(v);
    }

    const counts = {
      Parsimony: series.Parsimony.length,
      Dirichlet: series.Dirichlet.length,
      HSS: series.HSS.length,
    };

    const payload: ApiResp = { job_id: job, roots, root, series, counts, metric, layer };
    return NextResponse.json(payload);
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Internal error";
    return NextResponse.json({ error: message }, { status: 500 });
  }
}
