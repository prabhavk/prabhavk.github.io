// app/api/violin/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type MethodName = "Parsimony" | "Dirichlet" | "SSH";
type DBRow = { method: string; ll_final: number };

function normalizeMethod(s: string): MethodName | null {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "Parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "Dirichlet";
  if (m.includes("ssh")) return "SSH";
  return null;
}

type ApiResp = {
  job_id: string;
  root?: string;
  series: Record<MethodName, number[]>;
  counts: Record<MethodName, number>;
};

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job = (searchParams.get("job") ?? "").trim();
    const rootParam = (searchParams.get("root") ?? "").trim();
    const root = rootParam.length ? rootParam : undefined; // treat empty as undefined

    if (!job) {
      return NextResponse.json({ error: "Missing ?job=<job_id>" }, { status: 400 });
    }

    const params: string[] = [job];
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

    const payload: ApiResp = { job_id: job, root, series, counts };
    return NextResponse.json(payload);
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Internal error";
    return NextResponse.json({ error: message }, { status: 500 });
  }
}
