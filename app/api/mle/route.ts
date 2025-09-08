// app/api/mle/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type MleRow = {
  rep: string | number | null;
  // aliased names coming from SQL
  root_prob: string | number[] | null;
  trans_prob: string | number[][] | number[] | Record<string, unknown> | null;
  root: string | null;
  d_pi_1: number | string | null; d_pi_2: number | string | null; d_pi_3: number | string | null; d_pi_4: number | string | null;
  d_m_1:  number | string | null; d_m_2:  number | string | null; d_m_3:  number | string | null; d_m_4:  number | string | null;
};

type RepOnly = { rep: string | number | null };

/* ---------------- helpers ---------------- */

function normalizeMethod(s: string): "dirichlet" | "parsimony" | "hss" {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "dirichlet";
  if (m.includes("hss") || m.includes("hss")) return "hss";
  // default to dirichlet if empty/unknown
  return "dirichlet";
}

function tryParseJSON<T>(v: unknown): T | null {
  if (v == null) return null;
  if (typeof v === "string") { try { return JSON.parse(v) as T; } catch { return null; } }
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

/* ---------------- route ---------------- */

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

    if (!jobId) {
      return NextResponse.json({ error: "Missing job_id" }, { status: 400 });
    }

    const conn = db();

    // 1) list of reps for this job+method (numeric ordering even if stored as text)
    const repsSql = `
      SELECT DISTINCT br.rep
        FROM emtr_all_info AS br
       WHERE br.job_id = ? AND LOWER(br.method) = ?
       ORDER BY CAST(br.rep AS SIGNED)
    `;
    const { rows: repsRows } = await conn.execute<RepOnly>(repsSql, [jobId, method]);
    const reps: number[] = (repsRows || [])
      .map((r) => (r.rep == null ? NaN : Number(r.rep)))
      .filter((n) => Number.isFinite(n)) as number[];

    // 2) latest row for job+method (optionally for a specific rep)
    const rowSql = `
      SELECT
        br.rep,
        br.root_prob_final  AS root_prob,
        br.trans_prob_final AS trans_prob,
        br.root,
        j.d_pi_1, j.d_pi_2, j.d_pi_3, j.d_pi_4,
        j.d_m_1,  j.d_m_2,  j.d_m_3,  j.d_m_4
      FROM emtr_all_info AS br
      LEFT JOIN emtr_jobs AS j ON j.job_id = br.job_id
      WHERE br.job_id = ?
        AND LOWER(br.method) = ?
        ${repParam !== null ? "AND br.rep = ?" : ""}
      ORDER BY br.created_at DESC
      LIMIT 1
    `;
    const sqlArgs = repParam !== null
      ? [jobId, method, useNumericRepParam ? repParamNum : repParam]
      : [jobId, method];

    const { rows } = await conn.execute<MleRow>(rowSql, sqlArgs);

    if (!rows || rows.length === 0) {
      return NextResponse.json({
        job_id: jobId,
        method,
        rep: repParam !== null ? (useNumericRepParam ? repParamNum : repParam) : null,
        reps,
        root_prob: null,
        trans_prob: null,
        root: null,
        D_pi: null,
        D_M: null,
        // back-compat
        root_prob_final: null,
        trans_prob_final: null,
        trans: null,
      });
    }

    const row = rows[0];

    // Parse JSON/vector columns if they arrived as strings
    const root_prob = tryParseJSON<number[]>(row.root_prob);
    const transAny = tryParseJSON<unknown>(row.trans_prob);

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
      rep: repOut,
      reps,
      // new canonical names
      root_prob: root_prob ?? null,
      trans_prob: transAny ?? null,
      root: row.root ?? null,

      // back-compat fields
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
