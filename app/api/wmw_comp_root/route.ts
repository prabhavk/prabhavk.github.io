// app/api/wmw_comp_root/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type MethodName = "Parsimony" | "Dirichlet" | "HSS";

type Matrix = number[][];
type ApiResp = {
  job_id: string;
  method: MethodName;
  alpha: number;
  nodes: string[];
  sizes: Record<string, number>;
  p_value: Matrix;  // one-tailed p for H1: column > row
  z: Matrix;        // signed z (positive => medians suggest col>row)
  nA: Matrix;       // sample sizes (row node)
  nB: Matrix;       // sample sizes (col node)
  metric: "final" | "init"; // indicates which column was used
};

// ---- helpers ----
function normMethod(s: string | null): MethodName | null {
  if (!s) return null;
  const m = s.trim().toLowerCase();
  if (m.startsWith("par") || m.includes("pars")) return "Parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "Dirichlet";
  if (m.includes("hss")) return "HSS";
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

// Safe, lazy import of @bunchmark/stats.mwu
async function mwuZ(A: number[], B: number[]): Promise<number> {
  try {
    const mod = await import("@bunchmark/stats");
    const { z } = mod.mwu(A, B);
    return Number.isFinite(z) ? z : 0;
  } catch {
    return 0;
  }
}

// One-tailed test for H1: B > A (column > row)
async function oneTailedColGreater(
  A: number[],
  B: number[]
): Promise<{ p: number; z: number }> {
  if (A.length === 0 || B.length === 0) return { p: 1, z: 0 };
  const z_raw = await mwuZ(A, B);
  const mA = median(A), mB = median(B);

  const p_one = clamp01(1 - phi(Math.abs(z_raw)));
  if (!Number.isFinite(mA) || !Number.isFinite(mB)) return { p: 1, z: 0 };

  if (mB > mA) {
    return { p: p_one, z: Math.abs(z_raw) };
  } else {
    return { p: 1, z: -Math.abs(z_raw) };
  }
}

export async function GET(req: NextRequest) {
  const url = new URL(req.url);
  const job = (url.searchParams.get("job") || "").trim();
  const methodParam = normMethod(url.searchParams.get("method"));
  const alpha = Number(url.searchParams.get("alpha") ?? 0.05);
  const metricParam = (url.searchParams.get("metric") || "final").trim().toLowerCase();
  const metric: "final" | "init" = metricParam === "init" ? "init" : "final";

  if (!job) {
    return NextResponse.json({ ok: false, error: "Missing ?job=<job_id>" }, { status: 400 });
  }
  if (!methodParam) {
    return NextResponse.json(
      { ok: false, error: "Missing or invalid ?method=Parsimony|Dirichlet|HSS" },
      { status: 400 }
    );
  }
  if (!Number.isFinite(alpha) || alpha <= 0 || alpha >= 1) {
    return NextResponse.json({ ok: false, error: "Invalid alpha (0 < alpha < 1)" }, { status: 400 });
  }

  // Choose column based on metric
  const col = metric === "init" ? "ll_init" : "ll_final";

  const sql = `
    SELECT root, CAST(${col} AS DOUBLE) AS ll
    FROM emtr_init_final
    WHERE job_id = ?
      AND LOWER(method) = ?
      AND ${col} IS NOT NULL
    ORDER BY root, rep
  `;
  const params: (string | number)[] = [job, methodParam.toLowerCase()];

  try {
    const rows = await query<{ root: string; ll: number }>(sql, params);

    // Group by root
    const byNode = new Map<string, number[]>();
    for (const r of rows) {
      if (!r.root) continue;
      const v = Number(r.ll);
      if (!Number.isFinite(v)) continue;
      const arr = byNode.get(r.root) ?? [];
      arr.push(v);
      byNode.set(r.root, arr);
    }

    const nodes = Array.from(byNode.keys()).sort();
    const n = nodes.length;

    // Pre-allocate matrices
    const p_value: Matrix = Array.from({ length: n }, () => Array(n).fill(1));
    const z: Matrix = Array.from({ length: n }, () => Array(n).fill(0));
    const nA: Matrix = Array.from({ length: n }, () => Array(n).fill(0));
    const nB: Matrix = Array.from({ length: n }, () => Array(n).fill(0));

    const series = nodes.map((name) => byNode.get(name)!);

    const tasks: Array<Promise<void>> = [];
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        if (i === j) continue;
        const A = series[i];
        const B = series[j];
        nA[i][j] = A.length;
        nB[i][j] = B.length;

        tasks.push(
          (async () => {
            const { p, z: zval } = await oneTailedColGreater(A, B);
            p_value[i][j] = p;
            z[i][j] = zval;
          })()
        );
      }
    }
    await Promise.all(tasks);

    const sizes: Record<string, number> = {};
    nodes.forEach((name, idx) => (sizes[name] = series[idx].length));

    const payload: ApiResp = {
      job_id: job,
      method: methodParam,
      alpha,
      nodes,
      sizes,
      p_value,
      z,
      nA,
      nB,
      metric,
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
