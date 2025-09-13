// app/api/wmw_comp_method/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type MethodName = "Parsimony" | "Dirichlet" | "HSS";
type PairKey = "parsimony_vs_dirichlet" | "parsimony_vs_hss" | "dirichlet_vs_hss";

type DBRow = { method: string; val: unknown };

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
  // accept both legacy "ssh" and new "hss"
  if (m.includes("hss") || m.includes("ssh")) return "HSS";
  return null;
}

function toFiniteNumber(v: unknown): number | null {
  const n = Number(v);
  return Number.isFinite(n) ? n : null;
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

// safe MWU import without implicit any
async function safeMWU(A: number[], B: number[]): Promise<{ z: number }> {
  try {
    const mod: unknown = await import("@bunchmark/stats");
    const mwu = (mod as { mwu?: (a: number[], b: number[]) => { z: number } }).mwu;
    if (typeof mwu === "function") {
      const { z } = mwu(A, B);
      return { z: Number.isFinite(z) ? z : 0 };
    }
    return { z: 0 };
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
  // one-tailed p-value in the direction of the observed z
  const p_one = clamp01(1 - phi(Math.abs(z)));

  // choose winner by medians (stable even if z≈0)
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
  const job = (url.searchParams.get("job") || "").trim();
  const node = (url.searchParams.get("node") || "").trim() || undefined;
  const alpha = Number(url.searchParams.get("alpha") ?? 0.05);
  const metricParam = (url.searchParams.get("metric") || "final").trim().toLowerCase();
  const metric: "final" | "init" = metricParam === "init" ? "init" : "final"; // whitelist

  if (!job) {
    return NextResponse.json({ ok: false, error: "Missing ?job=<job_id>" }, { status: 400 });
  }
  if (!Number.isFinite(alpha) || alpha <= 0 || alpha >= 1) {
    return NextResponse.json({ ok: false, error: "Invalid alpha (0 < alpha < 1)" }, { status: 400 });
  }

  // whitelist the column to avoid SQL injection
  const col = metric === "init" ? "ll_init" : "ll_final";

  const nodeFilter = node ? "AND r.root = ?" : "";
  const sql = `
    SELECT LOWER(r.method) AS method, CAST(r.${col} AS DOUBLE) AS val
      FROM emtr_init_final r
      JOIN emtr_jobs j ON j.job_id = r.job_id
     WHERE r.job_id = ?
       ${nodeFilter}
       AND r.${col} IS NOT NULL
       AND j.status = 'completed'
     ORDER BY r.root, r.method, r.rep
  `;
  const params: (string | number)[] = node ? [job, node] : [job];

  try {
    const rows = await query<DBRow>(sql, params);

    const groups: Record<MethodName, number[]> = { Parsimony: [], Dirichlet: [], HSS: [] };
    for (const r of rows) {
      const m = normalizeMethod(r.method);
      const v = toFiniteNumber(r.val);
      if (!m || v === null) continue;
      groups[m].push(v);
    }

    const [pd, ps, ds] = await Promise.all([
      decideOneSided(groups.Parsimony, groups.Dirichlet, "Parsimony", "Dirichlet", alpha),
      decideOneSided(groups.Parsimony, groups.HSS,       "Parsimony", "HSS",       alpha),
      decideOneSided(groups.Dirichlet, groups.HSS,       "Dirichlet", "HSS",       alpha),
    ]);

    const payload: ApiResp = {
      job_id: job,
      node,
      alpha,
      sizes: {
        Parsimony: groups.Parsimony.length,
        Dirichlet: groups.Dirichlet.length,
        HSS: groups.HSS.length,
      },
      pairs: {
        parsimony_vs_dirichlet: pd,
        parsimony_vs_hss: ps,
        dirichlet_vs_hss: ds,
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
