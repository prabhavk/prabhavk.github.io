// app/api/mwu/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";
import { mwu } from "@bunchmark/stats";

export const runtime = "nodejs";

type MethodName = "Parsimony" | "Dirichlet" | "SSH";
type PairKey = "parsimony_vs_dirichlet" | "parsimony_vs_ssh" | "dirichlet_vs_ssh";

type DBRow = { method: string; ll_final: number };

type PairResult = {
  winner: MethodName | "none";
  p_value: number;   // one-sided p (winner's direction)
  z: number;         // z-statistic from MWU
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

// --- math helpers ---
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
function clamp01(v: number): number {
  if (!Number.isFinite(v)) return 1;
  if (v < 0) return 0;
  if (v > 1) return 1;
  return v;
}
function median(xs: number[]): number {
  const a = xs.slice().sort((x, y) => x - y);
  const n = a.length;
  if (n === 0) return NaN;
  const mid = Math.floor(n / 2);
  return n % 2 ? a[mid] : 0.5 * (a[mid - 1] + a[mid]);
}

// Decide winner with one-sided p = 1 - Phi(|z|); direction chosen from medians
function decideOneSided(
  A: number[], B: number[],
  nameA: MethodName, nameB: MethodName,
  alpha: number
): PairResult {
  if (A.length === 0 || B.length === 0) {
    return { winner: "none", p_value: 1, z: 0, nA: A.length, nB: B.length };
  }

  const { z } = mwu(A, B);
  const p_one = clamp01(1 - phi(Math.abs(z)));

  const mA = median(A);
  const mB = median(B);

  if (p_one < alpha) {
    const winner = mA > mB ? nameA : nameB;
    return { winner, p_value: p_one, z, nA: A.length, nB: B.length };
  }
  return { winner: "none", p_value: p_one, z, nA: A.length, nB: B.length };
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

    // If your node column isn't named `root`, change below.
    const params: string[] = node ? [job, node] : [job];
    const nodeFilter = node ? "AND root = ?" : "";

    const rows = await query<DBRow>(
      `
      SELECT method, ll_final
        FROM emtr_llchange
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

    const pairs: Record<PairKey, PairResult> = {
      parsimony_vs_dirichlet: decideOneSided(groups.Parsimony, groups.Dirichlet, "Parsimony", "Dirichlet", alpha),
      parsimony_vs_ssh:       decideOneSided(groups.Parsimony, groups.SSH,       "Parsimony", "SSH",       alpha),
      dirichlet_vs_ssh:       decideOneSided(groups.Dirichlet, groups.SSH,       "Dirichlet", "SSH",       alpha),
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
  } catch (err) {
    const message = err instanceof Error ? err.message : "Internal error";
    return NextResponse.json({ error: message }, { status: 500 });
  }
}
