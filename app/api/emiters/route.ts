// app/api/iters/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type MethodName = "Parsimony" | "Dirichlet" | "HSS";

type DBRow = {
  root: string;
  method: string;
  n: number;
  avg_iter: number;
  max_iter: number;
};

type Stats = { n: number; avg: number; max: number };
type ApiResp = {
  job_id: string;
  roots: string[];
  stats: Record<string, Record<MethodName, Stats>>; // root -> method -> stats
};

function normalizeMethod(s: string): MethodName | null {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "Parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "Dirichlet";
  if (m.includes("hss")) return "HSS";
  return null;
}

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job = searchParams.get("job") ?? "";

    if (!job) {
      return NextResponse.json({ error: "Missing ?job=<job_id>" }, { status: 400 });
    }

    // Aggregate in SQL
    const rows = await query<DBRow>(
      `
      SELECT root,
             method,
             COUNT(*)           AS n,
             AVG(iter)          AS avg_iter,
             MAX(iter)          AS max_iter
        FROM emtr_init_final
       WHERE job_id = ?
         AND iter IS NOT NULL
       GROUP BY root, method
      `,
      [job]
    );

    // Collect all roots and build stats map
    const rootsSet = new Set<string>();
    const stats: ApiResp["stats"] = {};

    for (const r of rows) {
      const root = String(r.root);
      const meth = normalizeMethod(r.method);
      if (!meth) continue;

      rootsSet.add(root);
      stats[root] ||= { Parsimony: { n: 0, avg: NaN, max: NaN }, Dirichlet: { n: 0, avg: NaN, max: NaN }, HSS: { n: 0, avg: NaN, max: NaN } };
      stats[root][meth] = {
        n: Number(r.n) || 0,
        avg: Number(r.avg_iter),
        max: Number(r.max_iter),
      };
    }

    // Sort roots like h_21..h_37 if possible
    const roots = Array.from(rootsSet).sort((a, b) => {
      const na = Number(a.replace(/^\D+/, "")); // extract numeric
      const nb = Number(b.replace(/^\D+/, ""));
      if (Number.isFinite(na) && Number.isFinite(nb)) return na - nb;
      return a.localeCompare(b);
    });

    const payload: ApiResp = { job_id: job, roots, stats };
    return NextResponse.json(payload);
  } catch (err: unknown) {
    const msg = err instanceof Error ? err.message : "Internal error";
    return NextResponse.json({ error: msg }, { status: 500 });
  }
}
