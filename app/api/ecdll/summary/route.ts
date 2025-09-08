// app/api/ecdll/summary/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type DbRow = {
  method: string;
  rep: number | null;
  root: string | null;
  root_prob_final: string | number[] | null;
  initial_ll: number | null; // aliased from ll_init
  final_ll: number | null;   // aliased from ll_final
  ecd_ll_per_iter: string | [number, number][] | Record<string, number> | null;
};

type Pt = [number, number];

/* ---------------- helpers ---------------- */

function normalizePointsUnknown(x: unknown): Pt[] {
  // Accept:
  // - [[iter, ll], ...]
  // - [{iter:..., ll:...}, ...]
  // - { "1": ll1, "2": ll2, ... }
  if (x == null) return [];
  const out: Pt[] = [];

  // Array case
  if (Array.isArray(x)) {
    for (const row of x) {
      if (Array.isArray(row) && row.length >= 2) {
        const it = Number(row[0]);
        const ll = Number(row[1]);
        if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
      } else if (row && typeof row === "object") {
        const r = row as Record<string, unknown>;
        const it = Number(r.iter);
        const ll = Number(r.ll);
        if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
      }
    }
    out.sort((a, b) => a[0] - b[0]);
    return out;
  }

  // Object map case: { "1": ll, "2": ll, ... }
  if (typeof x === "object") {
    for (const [k, v] of Object.entries(x as Record<string, unknown>)) {
      const it = Number(k);
      const ll = Number(v as number | string);
      if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
    }
    out.sort((a, b) => a[0] - b[0]);
    return out;
  }

  return [];
}

function normalizePoints(jsonish: string | [number, number][] | Record<string, number> | null): Pt[] {
  if (jsonish == null) return [];
  if (typeof jsonish === "string") {
    try {
      const parsed = JSON.parse(jsonish) as unknown;
      return normalizePointsUnknown(parsed);
    } catch {
      return [];
    }
  }
  return normalizePointsUnknown(jsonish);
}

function parseRoot(x: string | number[] | null): number[] | null {
  if (Array.isArray(x) && x.length === 4) {
    const nums = x.map((v) => Number(v));
    return nums.every((n) => Number.isFinite(n)) ? (nums as number[]) : null;
  }
  if (typeof x === "string") {
    try {
      const arr = JSON.parse(x) as unknown;
      if (Array.isArray(arr) && arr.length === 4) {
        const nums = arr.map((v) => Number(v));
        return nums.every((n) => Number.isFinite(n)) ? (nums as number[]) : null;
      }
    } catch {
      /* ignore */
    }
  }
  return null;
}

/* ---------------- route ---------------- */

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || "").trim();
    const repStr = searchParams.get("rep");
    const rep = repStr != null && repStr !== "" ? Number(repStr) : null;

    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }

    const conn = db();

    const baseSQL = `
      SELECT
        method,
        rep,
        root,
        root_prob_final,
        ll_init  AS initial_ll,
        ll_final AS final_ll,
        ecd_ll_per_iter
      FROM emtr_all_info
      WHERE job_id = ?
    `;

    const sql = rep !== null && Number.isFinite(rep)
      ? `${baseSQL} AND rep = ? ORDER BY created_at DESC`
      : `${baseSQL} ORDER BY created_at DESC`;

    const args = rep !== null && Number.isFinite(rep) ? [jobId, rep] : [jobId];
    const { rows } = await conn.execute<DbRow>(sql, args);

    // Keep latest per method only (thanks to ORDER BY created_at DESC)
    const latestByMethod = new Map<string, DbRow>();
    for (const r of rows) {
      if (!latestByMethod.has(r.method)) latestByMethod.set(r.method, r);
    }

    const methods = ["dirichlet", "parsimony", "hss"] as const;
    const out = methods.map((m) => {
      const row = latestByMethod.get(m);
      if (!row) {
        return {
          ok: true as const,
          method: m,                    
          root: null as number[] | null,
          initial_ll: null as number | null,
          final_ll: null as number | null,
          ecd_ll_first: null as number | null,
          ecd_ll_final: null as number | null,
          num_iterations: 0,
        };
      }

      const root = parseRoot(row.root_prob_final);
      const points = normalizePoints(row.ecd_ll_per_iter);
      const ecd_ll_first = points.length ? points[0][1] : null;
      const ecd_ll_final = points.length ? points[points.length - 1][1] : null;
      const num_iterations = points.length;

      return {
        ok: true as const,
        method: m,
        root,
        initial_ll: typeof row.initial_ll === "number" ? row.initial_ll : null,
        final_ll: typeof row.final_ll === "number" ? row.final_ll : null,
        ecd_ll_first,
        ecd_ll_final,
        num_iterations,
      };
    });

    return NextResponse.json({ ok: true, job_id: jobId, rep: rep ?? null, rows: out });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
