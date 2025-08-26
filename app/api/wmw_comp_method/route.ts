// app/api/wmw_comp_method/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

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

// --- helpers ---
function normalizeMethod(s: string): MethodName | null {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "Parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "Dirichlet";
  if (m.includes("ssh")) return "SSH";
  return null;
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
function phi(x: number): number {
  return 0.5 * (1 + erf(x / Math.SQRT2));
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

// safe MWU import
async function safeMWU(A: number[], B: number[]): Promise<{ z: number }> {
  try {
    const mod = await import("@bunchmark/stats");
    const { z } = mod.mwu(A, B);
    return { z: Number.isFinite(z) ? z : 0 };
  } catch {
    return { z: 0 };
  }
}

async function decideOneSided(
  A: number[], B: number[],
  nameA: MethodName, nameB: MethodName,
  alpha: number
): Promise<PairResult> {
  if (A.length === 0 || B.length === 0) {
    return { winner: "none", p_value: 1, z: 0, nA: A.length, nB: B.length };
  }

  const { z } = await safeMWU(A, B);
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
  const url = new URL(req.url);
  const job = url.searchParams.get("job") || "";
  const node = url.searchParams.get("node") || undefined;
  const alpha = Number(url.searchParams.get("alpha") ?? 0.05);

  if (!job) {
    return NextResponse.json({ ok: false, error: "Missing ?job=<job_id>" }, { status: 400 });
  }
  if (!Number.isFinite(alpha) || alpha <= 0 || alpha >= 1) {
    return NextResponse.json({ ok: false, error: "Invalid alpha (0 < alpha < 1)" }, { status: 400 });
  }

  const nodeFilter = node ? "AND root = ?" : "";
  const sql = `
    SELECT LOWER(method) AS method, CAST(ll_final AS DOUBLE) AS ll_final
    FROM emtr_init_final
    WHERE job_id = ?
      ${nodeFilter}
      AND ll_final IS NOT NULL
    ORDER BY root, method, rep
  `;
  const params: (string | number)[] = node ? [job, node] : [job];

  try {
    const rows = await query<DBRow>(sql, params);

    const groups: Record<MethodName, number[]> = { Parsimony: [], Dirichlet: [], SSH: [] };
    for (const r of rows) {
      const m = normalizeMethod(r.method);
      if (!m) continue;
      const v = Number(r.ll_final);
      if (Number.isFinite(v)) groups[m].push(v);
    }

    const [pd, ps, ds] = await Promise.all([
      decideOneSided(groups.Parsimony, groups.Dirichlet, "Parsimony", "Dirichlet", alpha),
      decideOneSided(groups.Parsimony, groups.SSH,       "Parsimony", "SSH",       alpha),
      decideOneSided(groups.Dirichlet, groups.SSH,       "Dirichlet", "SSH",       alpha),
    ]);

    const payload: ApiResp = {
      job_id: job,
      node,
      alpha,
      sizes: {
        Parsimony: groups.Parsimony.length,
        Dirichlet: groups.Dirichlet.length,
        SSH: groups.SSH.length,
      },
      pairs: {
        parsimony_vs_dirichlet: pd,
        parsimony_vs_ssh: ps,
        dirichlet_vs_ssh: ds,
      },
    };

    return NextResponse.json(payload, { status: 200 });
  } catch (err) {
    const message = err instanceof Error ? err.message : String(err);
    return NextResponse.json(
      { ok: false, error: "Query failed", detail: message, sql: sql.trim(), params },
      { status: 500 }
    );
  }
}
