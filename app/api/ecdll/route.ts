// app/api/ecdll/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Row = {
  ecd_ll_per_iter: unknown;
  root: string | null;
  root_prob_final: unknown;
  ll_init: number | string | null;
  ll_final: number | string | null;
  rep: number | string | null;
};

type Pt = [number, number];

/* ---------------- helpers ---------------- */

function normalizePointsUnknown(x: unknown): Pt[] {
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
        const itRaw = "iter" in obj ? (obj as Record<string, unknown>).iter : undefined;
        const llRaw = "ll"   in obj ? (obj as Record<string, unknown>).ll   : undefined;

        const it = Number(itRaw);
        const ll = Number(llRaw);

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
        const nums = v.map((x: unknown) => Number(x));
        return nums.every((n) => Number.isFinite(n)) ? (nums as number[]) : null;
      }
    } catch { /* ignore */ }
  }
  return null;
}

const toNum = (x: unknown): number | null => {
  const n = Number(x);
  return Number.isFinite(n) ? n : null;
};

/* ---------------- route ---------------- */

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || searchParams.get("job") || "").trim();

    // method provided by the client; we normalize to lowercase for comparisons
    const methodRaw = (searchParams.get("method") || "").trim();
    const methodLc = methodRaw.toLowerCase();

    // optional rep and root filters
    const repStr = searchParams.get("rep");
    const rep = repStr != null && repStr !== "" ? Number(repStr) : null;

    const rootParam = (searchParams.get("root") || "").trim(); // <-- NEW: root filter

    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }
    if (!methodRaw) {
      return NextResponse.json({ ok: false, error: "Missing method" }, { status: 400 });
    }

    const conn = db();

    // Base WHERE (we’ll append rep/root conditions below as needed)
    const where = `WHERE job_id = ? AND LOWER(method) = ?`;
    const baseParams: Array<string | number> = [jobId, methodLc];

    // Helper to run the “latest row” query with optional extra WHERE + params
    const runLatest = async (extraWhere: string, extraParams: Array<string | number>) => {
      const sql = `
        SELECT
          ecd_ll_per_iter,
          root,
          root_prob_final,
          ll_init,
          ll_final,
          rep
        FROM emtr_all_info
        ${where} ${extraWhere}
        ORDER BY created_at DESC
        LIMIT 1
      `;
      const { rows } = await conn.execute<Row>(sql, [...baseParams, ...extraParams]);
      return rows ?? [];
    };

    let rows: Row[] = [];

    // Priority 1: exact (rep, root) match if both provided
    if (rep !== null && Number.isFinite(rep) && rootParam) {
      rows = await runLatest(`AND rep = ? AND root = ?`, [rep, rootParam]);
      // Fallback if not found: same rep, ignore root
      if (!rows.length) rows = await runLatest(`AND rep = ?`, [rep]);
      // Fallback if still not found: ignore rep, keep root
      if (!rows.length) rows = await runLatest(`AND root = ?`, [rootParam]);
    }
    // Priority 2: only rep provided
    else if (rep !== null && Number.isFinite(rep)) {
      rows = await runLatest(`AND rep = ?`, [rep]);
      // Fallback: latest for this method (and root, if provided)
      if (!rows.length && rootParam) rows = await runLatest(`AND root = ?`, [rootParam]);
    }
    // Priority 3: only root provided
    else if (rootParam) {
      rows = await runLatest(`AND root = ?`, [rootParam]);
    }

    // Final fallback: latest regardless of rep/root for this method
    if (!rows.length) {
      rows = await runLatest(``, []);
    }

    if (!rows || rows.length === 0) {
      return NextResponse.json({
        ok: true,
        job_id: jobId,
        method: methodLc,
        rep: rep ?? null,
        points: [],
        summary: null,
      });
    }

    const row = rows[0];

    const points = parsePoints(row.ecd_ll_per_iter);
    const ecd_first = points.length ? points[0][1] : null;
    const ecd_final = points.length ? points[points.length - 1][1] : null;

    // Coerce numerics in case driver hands back strings
    const repOut = row.rep == null ? null : toNum(row.rep);
    const llInit = row.ll_init == null ? null : toNum(row.ll_init);
    const llFinal = row.ll_final == null ? null : toNum(row.ll_final);

    return NextResponse.json({
      ok: true,
      job_id: jobId,
      method: methodLc,
      rep: repOut ?? rep ?? null,
      points,
      summary: {
        root: row.root ?? null,
        num_iterations: points.length,
        ll_init: llInit,
        ecd_ll_first: ecd_first,
        ecd_ll_final: ecd_final,
        ll_final: llFinal,
        root_prob_final: parseRootProb(row.root_prob_final),
      },
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
