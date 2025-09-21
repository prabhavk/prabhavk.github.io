// app/api/mle/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type MleRow = {
  rep: string | number | null;
  root: string | null;
  // from emtr_trifle_layers
  root_prob: string | number[] | null;
  trans_prob: string | number[][] | number[] | Record<string, unknown> | null;
  // hyperparams from emtr_jobs
  d_pi_1: number | string | null; d_pi_2: number | string | null; d_pi_3: number | string | null; d_pi_4: number | string | null;
  d_m_1:  number | string | null; d_m_2:  number | string | null; d_m_3:  number | string | null; d_m_4:  number | string | null;
};

type RepOnly = { rep: string | number | null };

function normalizeMethod(s: string): "dirichlet" | "parsimony" | "hss" {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "dirichlet";
  if (m.includes("hss")) return "hss";
  return "dirichlet";
}

function tryParseJSON<T>(v: unknown): T | null {
  if (v == null) return null;
  if (typeof v === "string") {
    try {
      return JSON.parse(v) as T;
    } catch {
      return null;
    }
  }
  return v as T;
}

function numOrNaN(v: unknown): number {
  const n = typeof v === "string" ? Number(v.trim()) : Number(v);
  return Number.isFinite(n) ? n : NaN;
}

function toNumArrayOrNull(a: Array<number | string | null | undefined>): number[] | null {
  const nums = a.map(numOrNaN);
  return nums.every((x) => Number.isFinite(x) && x > 0) ? (nums as number[]) : null;
}

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || "").trim();
    const methodRaw = (searchParams.get("method") || "dirichlet").trim();
    const method = normalizeMethod(methodRaw);

    const repParamStr = (searchParams.get("rep") || "").trim();
    const repParam = repParamStr !== "" ? repParamStr : null;
    const repParamNum = repParam !== null ? Number(repParam) : null;
    const useNumericRepParam = repParamNum !== null && Number.isFinite(repParamNum);

    const layerStr = (searchParams.get("layer") || "0").trim();
    const layerNum = Number(layerStr);
    const layer = layerNum === 0 || layerNum === 1 || layerNum === 2 ? layerNum : 0;

    const rootFilter = (searchParams.get("root") || "").trim();
    const hasRootFilter = rootFilter.length > 0;

    if (!jobId) {
      return NextResponse.json({ error: "Missing job_id" }, { status: 400 });
    }

    const conn = db();

    // 1) Distinct reps for this job+method from emtr_trifle_layers
    const repsSql = `
      SELECT DISTINCT tl.rep
        FROM emtr_trifle_layers AS tl
       WHERE tl.job_id = ? AND LOWER(tl.method) = ?
       ORDER BY CAST(tl.rep AS SIGNED)
    `;
    const { rows: repsRows } = await conn.execute<RepOnly>(repsSql, [jobId, method]);
    const reps: number[] = (repsRows || [])
      .map((r) => (r.rep == null ? NaN : Number(r.rep)))
      .filter((n) => Number.isFinite(n)) as number[];

    // 2) Latest matching row for job+method+layer (optional rep/root)
    const rowSql = `
      SELECT
        tl.rep,
        tl.root,
        tl.root_prob_final  AS root_prob,
        tl.trans_prob_final AS trans_prob,
        j.d_pi_1, j.d_pi_2, j.d_pi_3, j.d_pi_4,
        j.d_m_1,  j.d_m_2,  j.d_m_3,  j.d_m_4
      FROM emtr_trifle_layers AS tl
      LEFT JOIN emtr_jobs AS j ON j.job_id = tl.job_id
      WHERE tl.job_id = ?
        AND LOWER(tl.method) = ?
        AND tl.layer = ?
        ${repParam !== null ? "AND tl.rep = ?" : ""}
        ${hasRootFilter ? "AND tl.root = ?" : ""}
      ORDER BY tl.created_at DESC
      LIMIT 1
    `;
    const sqlArgs: Array<string | number> = [jobId, method, layer];
    if (repParam !== null) sqlArgs.push(useNumericRepParam ? (repParamNum as number) : repParam);
    if (hasRootFilter) sqlArgs.push(rootFilter);

    const { rows } = await conn.execute<MleRow>(rowSql, sqlArgs);

    if (!rows || rows.length === 0) {
      return NextResponse.json({
        job_id: jobId,
        method,
        layer,
        rep: repParam !== null ? (useNumericRepParam ? repParamNum : repParam) : null,
        reps,
        root_prob: null,
        trans_prob: null,
        root: hasRootFilter ? rootFilter : null,
        D_pi: null,
        D_M: null,
        // back-compat
        root_prob_final: null,
        trans_prob_final: null,
        trans: null,
      });
    }

    const row = rows[0];

    const root_prob = tryParseJSON<number[]>(row.root_prob);
    const transAny  = tryParseJSON<unknown>(row.trans_prob);

    // Dirichlet hyperparams (accept strings or numbers)
    const D_pi = toNumArrayOrNull([row.d_pi_1, row.d_pi_2, row.d_pi_3, row.d_pi_4]);
    const D_M  = toNumArrayOrNull([row.d_m_1,  row.d_m_2,  row.d_m_3,  row.d_m_4]);

    const repOut =
      row.rep == null
        ? null
        : (Number.isFinite(Number(row.rep)) ? Number(row.rep) : String(row.rep));

    return NextResponse.json({
      job_id: jobId,
      method,
      layer,
      rep: repOut,
      reps,
      // canonical
      root_prob: root_prob ?? null,
      trans_prob: transAny ?? null,
      root: hasRootFilter ? rootFilter : (row.root ?? null),

      // back-compat
      root_prob_final: root_prob ?? null,
      trans_prob_final: transAny ?? null,
      trans: (transAny as unknown) ?? null,

      D_pi,
      D_M,
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ error: msg }, { status: 500 });
  }
}
