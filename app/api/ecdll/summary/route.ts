// app/api/ecdll/summary/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type DbRow = {
  method: string;
  root_prob_final: string | number[] | null;
  initial_ll: number | null;
  final_ll: number | null;
  ecd_ll_per_iter: string | [number, number][] | null;
};

type Pt = [number, number];

function normalizePoints(x: unknown): Pt[] {
  if (!Array.isArray(x)) return [];
  const out: Pt[] = [];
  for (const row of x) {
    if (Array.isArray(row) && row.length >= 2) {
      const it = Number(row[0]);
      const ll = Number(row[1]);
      if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
    } else if (row && typeof row === "object") {
      const it = Number((row as Record<string, unknown>).iter);
      const ll = Number((row as Record<string, unknown>).ll);
      if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
    }
  }
  out.sort((a, b) => a[0] - b[0]);
  return out;
}

function parseRoot(x: string | number[] | null): number[] | null {
  if (Array.isArray(x) && x.length === 4) return x.map(Number);
  if (typeof x === "string") {
    try {
      const arr = JSON.parse(x) as unknown;
      if (Array.isArray(arr) && arr.length === 4) return arr.map((v) => Number(v));
    } catch {
      /* ignore */
    }
  }
  return null;
}

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || "").trim();
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }

    const sql = `
      SELECT method, root_name, ll_init, ll_final, 
      FROM emtr_best_rep
      WHERE job_id = ?
      ORDER BY created_at DESC
    `;

    const conn = db();
    const { rows } = await conn.execute<DbRow>(sql, [jobId]);

    // Keep latest per method only
    const latestByMethod = new Map<string, DbRow>();
    for (const r of rows) {
      if (!latestByMethod.has(r.method)) latestByMethod.set(r.method, r);
    }

    const methods = ["dirichlet", "parsimony", "ssh"] as const;
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
      const points =
        typeof row.ecd_ll_per_iter === "string"
          ? normalizePoints(JSON.parse(row.ecd_ll_per_iter) as unknown)
          : normalizePoints(row.ecd_ll_per_iter);

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

    return NextResponse.json({ ok: true, job_id: jobId, rows: out });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
