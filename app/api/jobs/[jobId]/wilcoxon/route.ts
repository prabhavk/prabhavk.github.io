// app/api/jobs/[jobId]/wilcoxon/route.ts
import { NextResponse } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Method = "parsimony" | "dirichlet" | "ssh";
type Row = {
  root: string;
  method: Method;
  rep: number;
  ll_final: number;
};

/* -------------------- stats (no deps) -------------------- */

// average ranks with ties; returns ranks and sizes of tie groups
function rankWithTies(values: number[]): { ranks: number[]; tieGroups: number[] } {
  const n = values.length;
  const idx = values.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
  const ranks = new Array<number>(n);
  const tieGroups: number[] = [];
  let i = 0;
  while (i < n) {
    let j = i + 1;
    while (j < n && idx[j].v === idx[i].v) j++;
    const size = j - i;
    const avgRank = (i + 1 + j) / 2;
    for (let k = i; k < j; k++) ranks[idx[k].i] = avgRank;
    if (size > 1) tieGroups.push(size);
    i = j;
  }
  return { ranks, tieGroups };
}

function erf(x: number) {
  const sign = Math.sign(x);
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
  const t = 1 / (1 + p * Math.abs(x));
  const y = 1 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * Math.exp(-x * x);
  return sign * y;
}
const normalCdf = (z: number) => 0.5 * (1 + erf(z / Math.SQRT2));

/** Mann–Whitney U with tie-corrected normal approximation; two-sided p */
function mannWhitney(xIn: number[], yIn: number[]) {
  const x = xIn.filter(Number.isFinite);
  const y = yIn.filter(Number.isFinite);
  const n1 = x.length, n2 = y.length;
  if (n1 < 1 || n2 < 1) return { U: NaN, p: NaN, n1, n2 };

  const all = x.concat(y);
  const N = n1 + n2;
  const { ranks, tieGroups } = rankWithTies(all);

  // sum of ranks for group X (first n1)
  let R1 = 0;
  for (let i = 0; i < n1; i++) R1 += ranks[i];

  const U1 = R1 - (n1 * (n1 + 1)) / 2;
  const U2 = n1 * n2 - U1;
  const U = Math.min(U1, U2);

  // variance with tie correction
  let tieSum = 0;
  for (const t of tieGroups) tieSum += t * t * t - t;
  const varU = (n1 * n2 / 12) * ((N + 1) - (tieSum / (N * (N - 1))));
  const sigma = Math.sqrt(Math.max(varU, 0));
  const mu = (n1 * n2) / 2;

  // continuity correction
  const z = sigma > 0 ? (U - mu + 0.5 * Math.sign(mu - U)) / sigma : 0;
  const p = 2 * (1 - normalCdf(Math.abs(z)));

  return { U, p, n1, n2 };
}

/* -------------------- context helpers (no 'any') -------------------- */

function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isPromise<T = unknown>(x: unknown): x is Promise<T> {
  return isObject(x) && "then" in x && typeof (x as { then?: unknown }).then === "function";
}
async function getJobIdFromContext(ctx: unknown): Promise<string | null> {
  if (!isObject(ctx) || !("params" in ctx)) return null;
  const p = (ctx as { params: unknown }).params;

  if (isPromise<Record<string, unknown>>(p)) {
    const resolved = await p;
    return isObject(resolved) && typeof resolved.jobId === "string" ? resolved.jobId : null;
  }
  return isObject(p) && typeof p.jobId === "string" ? p.jobId : null;
}

/* -------------------- route handler -------------------- */

export async function GET(_req: Request, ctx: unknown) {
  try {
    const jobId = await getJobIdFromContext(ctx);
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }

    const conn = db();
    const res = await conn.execute(
      `
      SELECT root, LOWER(method) AS method, rep, ll_final
      FROM emtr_init_final
      WHERE job_id = ?
        AND LOWER(method) IN ('parsimony','dirichlet','ssh')
      ORDER BY root, method, rep
      `,
      [jobId]
    );
    const rows = res.rows as Row[];

    const byRoot = new Map<string, Record<Method, number[]>>();
    for (const r of rows) {
      if (!byRoot.has(r.root)) byRoot.set(r.root, { parsimony: [], dirichlet: [], ssh: [] });
      byRoot.get(r.root)![r.method].push(r.ll_final);
    }

    const results = Array.from(byRoot, ([root, s]) => ({
      root,
      parsimony_vs_dirichlet: mannWhitney(s.parsimony, s.dirichlet),
      parsimony_vs_ssh:       mannWhitney(s.parsimony, s.ssh),
      dirichlet_vs_ssh:       mannWhitney(s.dirichlet, s.ssh),
    }));

    return NextResponse.json({ ok: true, jobId, results });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
