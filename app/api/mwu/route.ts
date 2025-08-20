// app/api/mwu/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";
import { mwu } from "@bunchmark/stats";

export const runtime = "nodejs";

type MethodName = "Parsimony" | "Dirichlet" | "SSH";
type PairKey = "parsimony_vs_dirichlet" | "parsimony_vs_ssh" | "dirichlet_vs_ssh";

type DBRow = { method: string; ll_final: number }; // adjust names below if your schema differs

type PairResult = {
  winner: MethodName | "none";
  p_value: number;   // one-sided p for the winning direction
  nA: number;
  nB: number;
};

type ApiResp = {
  job_id: string;
  node?: string;
  alpha: number;
  sizes: Record<MethodName, number>;
  pairs: Record<PairKey, PairResult>;
};

function normalizeMethod(s: string): MethodName | null {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "Parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "Dirichlet";
  if (m.includes("ssh")) return "SSH";
  return null;
}

// Normal CDF via erf — for one-sided p from z
function phi(x: number): number {
  return 0.5 * (1 + erf(x / Math.SQRT2));
}
function erf(x: number): number {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
  const sign = x < 0 ? -1 : 1;
  const ax = Math.abs(x);
  const t = 1 / (1 + p * ax);
  const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-ax * ax);
  return sign * y;
}

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job = searchParams.get("job");
    const node = searchParams.get("node") ?? undefined; // e.g., "h_21"
    const alpha = Number(searchParams.get("alpha") ?? 0.05);

    if (!job) return NextResponse.json({ error: "Missing ?job=<job_id>" }, { status: 400 });
    if (!Number.isFinite(alpha) || alpha <= 0 || alpha >= 1) {
      return NextResponse.json({ error: "Invalid alpha (0 < alpha < 1)" }, { status: 400 });
    }

    // NOTE: If your node column isn't named `root`, change `root = ?` to the correct column.
    const params: Array<string> = node ? [job, node] : [job];
    const nodeFilter = node ? "AND root = ?" : ""; // <-- rename `root` if needed

    const rows = await query<DBRow>(
      `
      SELECT method, ll_final
        FROM emtr_rows
       WHERE job_id = ?
         ${nodeFilter}
         AND ll_final IS NOT NULL
      `,
      params
    );

    const groups: Record<MethodName, number[]> = { Parsimony: [], Dirichlet: [], SSH: [] };
    for (const r of rows) {
      const m = normalizeMethod(r.method);
      if (!m) continue;
      const v = Number(r.ll_final);
      if (Number.isFinite(v)) groups[m].push(v);
    }

    const decide = (nameA: MethodName, nameB: MethodName): PairResult => {
      const A = groups[nameA], B = groups[nameB];
      if (A.length === 0 || B.length === 0) {
        return { winner: "none", p_value: 1, nA: A.length, nB: B.length };
      }
      const { z } = mwu(A, B); // from @bunchmark/stats
      // One-sided p-values for the two directions
      const pAgtB = 1 - phi(z); // H1: A > B
      const pBgtA = phi(z);     // H1: B > A

      if (pAgtB < alpha) return { winner: nameA, p_value: pAgtB, nA: A.length, nB: B.length };
      if (pBgtA < alpha) return { winner: nameB, p_value: pBgtA, nA: A.length, nB: B.length };
      return { winner: "none", p_value: Math.min(pAgtB, pBgtA), nA: A.length, nB: B.length };
    };

    const pairs: ApiResp["pairs"] = {
      parsimony_vs_dirichlet: decide("Parsimony", "Dirichlet"),
      parsimony_vs_ssh:       decide("Parsimony", "SSH"),
      dirichlet_vs_ssh:       decide("Dirichlet", "SSH"),
    };

    const payload: ApiResp = {
      job_id: job,
      node,
      alpha,
      sizes: {
        Parsimony: groups.Parsimony.length,
        Dirichlet: groups.Dirichlet.length,
        SSH: groups.SSH.length,
      },
      pairs,
    };

    return NextResponse.json(payload);
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Internal error";
    return NextResponse.json({ error: message }, { status: 500 });
  }
}
