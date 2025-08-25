// app/api/ecdll/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Row = {
  ecd_ll_per_iter: unknown;               // <- allow string | array | object
  root_name: string | null;
  root_prob_final: unknown;               // <- allow string | array
  ll_init: number | null;
  ll_final: number | null;
  rep: number | null;
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

  if (Array.isArray(x)) {
    for (const r of x) {
      if (Array.isArray(r) && r.length >= 2) {
        const it = Number(r[0]);
        const ll = Number(r[1]);
        if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
      } else if (r && typeof r === "object") {
        const obj = r as Record<string, unknown>;
        const it = Number(obj.iter);
        const ll = Number(obj.ll);
        if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
      }
    }
    out.sort((a, b) => a[0] - b[0]);
    return out;
  }

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

function parsePoints(jsonish: unknown): Pt[] {
  if (jsonish == null) return [];
  if (typeof jsonish === "string") {
    try {
      return normalizePointsUnknown(JSON.parse(jsonish) as unknown);
    } catch {
      return [];
    }
  }
  return normalizePointsUnknown(jsonish);
}

function parseRootProb(jsonish: unknown): number[] | null {
  if (Array.isArray(jsonish) && jsonish.length === 4) {
    const nums = jsonish.map((x) => Number(x));
    return nums.every((n) => Number.isFinite(n)) ? (nums as number[]) : null;
  }
  if (typeof jsonish === "string") {
    try {
      const v = JSON.parse(jsonish);
      if (Array.isArray(v) && v.length === 4) {
        const nums = v.map((x) => Number(x));
        return nums.every((n) => Number.isFinite(n)) ? (nums as number[]) : null;
      }
    } catch {/* ignore */}
  }
  return null;
}

/* ---------------- route ---------------- */

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || "").trim();
    const method = (searchParams.get("method") || "").trim().toLowerCase();
    const repStr = searchParams.get("rep");
    const rep = repStr != null && repStr !== "" ? Number(repStr) : null;

    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }
    if (!method) {
      return NextResponse.json({ ok: false, error: "Missing method" }, { status: 400 });
    }

    const baseSQL = `
      SELECT
        ecd_ll_per_iter,
        root_name,
        root_prob_final,
        ll_init,
        ll_final,
        rep
      FROM emtr_all_info
      WHERE job_id = ? AND method = ?
    `;

    const conn = db();

    let rows: Row[] = [];
    if (rep !== null && Number.isFinite(rep)) {
      // Try specific rep first
      const sql = `${baseSQL} AND rep = ? ORDER BY created_at DESC LIMIT 1`;
      rows = (await conn.execute<Row>(sql, [jobId, method, rep])).rows ?? [];
      // Fallback to latest if no exact rep row exists
      if (!rows.length) {
        rows = (await conn.execute<Row>(`${baseSQL} ORDER BY created_at DESC LIMIT 1`, [jobId, method])).rows ?? [];
      }
    } else {
      rows = (await conn.execute<Row>(`${baseSQL} ORDER BY created_at DESC LIMIT 1`, [jobId, method])).rows ?? [];
    }

    if (!rows || rows.length === 0) {
      return NextResponse.json({
        ok: true,
        job_id: jobId,
        method,
        rep: rep ?? null,
        points: [],
        summary: null,
      });
    }

    const row = rows[0];

    const points = parsePoints(row.ecd_ll_per_iter);
    const ecd_first = points.length ? points[0][1] : null;
    const ecd_final = points.length ? points[points.length - 1][1] : null;

    return NextResponse.json({
      ok: true,
      job_id: jobId,
      method,
      rep: row.rep ?? rep ?? null,
      points,
      summary: {
        root: row.root_name ?? null,
        num_iterations: points.length,
        ll_init: row.ll_init,
        ecd_ll_first: ecd_first,
        ecd_ll_final: ecd_final,
        ll_final: row.ll_final,
        root_prob_final: parseRootProb(row.root_prob_final),
      },
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
