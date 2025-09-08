// app/api/allinfo/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

/* -------------------- Types -------------------- */

type EMStruct = {
  method?: string;
  rep?: number | null;
  iter?: number | null;
  ecd_ll_per_iter?: Record<string, number> | number[] | null;
  ll_init?: number | null;
  ll_final?: number | null;
  root?: string | null;
  root_prob_init?: number[] | null;
  root_prob_final?: number[] | null;
  trans_prob_init?: Record<string, unknown> | null;
  trans_prob_final?: Record<string, unknown> | null;
};

type BestRow = {
  method: string | null;
  rep: number | null;
  iter: number | null;
  root: string | null;
  ll_init: number | null;
  ll_final: number | null;
  root_prob_init: string | null;
  root_prob_final: string | null;
  trans_prob_init: string | null;
  trans_prob_final: string | null;
  ecd_ll_per_iter: string | null;
  raw_json: string | null;
};

type JobParamRow = {
  d_pi_1: number | null; d_pi_2: number | null; d_pi_3: number | null; d_pi_4: number | null;
  d_m_1:  number | null; d_m_2:  number | null; d_m_3:  number | null; d_m_4:  number | null;
};

/* -------------------- Helpers -------------------- */

function asObj(x: unknown): Record<string, unknown> | null {
  return x && typeof x === "object" ? (x as Record<string, unknown>) : null;
}

function toStrOrNull(v: unknown): string | null {
  return v == null ? null : String(v);
}

function toNumOrNull(v: unknown): number | null {
  if (v == null) return null;
  const n = Number(v);
  return Number.isFinite(n) ? n : null;
}

function coerceBestRow(u: unknown): BestRow {
  const o = (u && typeof u === "object" ? (u as Record<string, unknown>) : {}) as Record<string, unknown>;
  return {
    method: toStrOrNull(o.method),
    rep: toNumOrNull(o.rep),
    iter: toNumOrNull(o.iter),
    root: toStrOrNull(o.root),
    ll_init: toNumOrNull(o.ll_init),
    ll_final: toNumOrNull(o.ll_final),
    root_prob_init: toStrOrNull(o.root_prob_init),
    root_prob_final: toStrOrNull(o.root_prob_final),
    trans_prob_init: toStrOrNull(o.trans_prob_init),
    trans_prob_final: toStrOrNull(o.trans_prob_final),
    ecd_ll_per_iter: toStrOrNull(o.ecd_ll_per_iter),
    raw_json: toStrOrNull(o.raw_json),
  };
}

function coerceJobParamRow(u: unknown): JobParamRow {
  const o = (u && typeof u === "object" ? (u as Record<string, unknown>) : {}) as Record<string, unknown>;
  return {
    d_pi_1: toNumOrNull(o.d_pi_1), d_pi_2: toNumOrNull(o.d_pi_2), d_pi_3: toNumOrNull(o.d_pi_3), d_pi_4: toNumOrNull(o.d_pi_4),
    d_m_1:  toNumOrNull(o.d_m_1),  d_m_2:  toNumOrNull(o.d_m_2),  d_m_3:  toNumOrNull(o.d_m_3),  d_m_4:  toNumOrNull(o.d_m_4),
  };
}

type ExecResult = { rows?: unknown[] };

/** Safely extract typed rows from a PlanetScale execute() result */
function rowsFromResult<T>(rs: ExecResult, map: (u: unknown) => T): T[] {
  const base = rs && typeof rs === "object" ? rs : null;
  if (!base || !Array.isArray(base.rows)) return [];
  return (base.rows as unknown[]).map(map);
}

function parseMaybeJson<T>(v: unknown): T | null {
  if (v == null) return null;
  if (typeof v === "string") {
    try { return JSON.parse(v) as T; } catch { return null; }
  }
  return (v as T) ?? null;
}

/** Accept worker payloads of different shapes and normalize to {pars?, dirichlet?, hss?}. */
function normalizeIncoming(obj: Record<string, unknown>): {
  pars?: EMStruct;
  dirichlet?: EMStruct;
  hss?: EMStruct;
} {
  const out: { pars?: EMStruct; dirichlet?: EMStruct; hss?: EMStruct } = {};

  const getStruct = (k: string): EMStruct | undefined => {
    const v = asObj(obj[k]);
    return v ? (v as unknown as EMStruct) : undefined;
  };

  // Common keyed variants
  const pars1 = getStruct("pars") ?? getStruct("parsimony");
  const dir1  = getStruct("dirichlet");
  const hss1  = getStruct("hss");

  if (pars1) out.pars = pars1;
  if (dir1)  out.dirichlet = dir1;
  if (hss1)  out.hss = hss1;

  // Single-method struct at top-level (no keys)
  if (!out.pars && !out.dirichlet && !out.hss) {
    const maybeSingle = obj as Partial<EMStruct>;
    const m = typeof maybeSingle.method === "string" ? maybeSingle.method.toLowerCase() : "";
    if (m.includes("pars")) out.pars = maybeSingle as EMStruct;
    else if (m.includes("dir")) out.dirichlet = maybeSingle as EMStruct;
    else if (m.includes("hss")) out.hss = maybeSingle as EMStruct;
  }

  return out;
}

/* -------------------- DB upsert -------------------- */

type MethodKey = "parsimony" | "dirichlet" | "hss";

function computeNumIter(ecd: EMStruct["ecd_ll_per_iter"]): number | null {
  if (!ecd) return null;
  if (Array.isArray(ecd)) return ecd.length;
  if (typeof ecd === "object") return Object.keys(ecd).length;
  return null;
}

async function upsertOne(
  conn: ReturnType<typeof db>,
  job_id: string,
  methodKey: MethodKey,
  em: EMStruct
) {
  // If your table’s PK is (job_id, method, rep) you’ll get one row per rep.
  // If your table’s PK is (job_id, method) only, the latest write overwrites the previous rep.
  const rep = Number.isFinite(em.rep as number) ? (em.rep as number) : 1;

  const root        = typeof em.root === "string" ? em.root : null;
  const ll_init          = Number.isFinite(em.ll_init as number) ? (em.ll_init as number) : null;
  const ll_final         = Number.isFinite(em.ll_final as number) ? (em.ll_final as number) : null;
  const iter         = Number.isFinite(em.iter as number) ? (em.iter as number) : computeNumIter(em.ecd_ll_per_iter);

  const root_prob_init   = em.root_prob_init ?? null;
  const root_prob_final  = em.root_prob_final ?? null;
  const trans_prob_init  = em.trans_prob_init ?? null;
  const trans_prob_final = em.trans_prob_final ?? null;
  const ecd_ll_per_iter  = em.ecd_ll_per_iter ?? null;

  const raw_json: unknown = em; // keep full struct

  const sql = `
    INSERT INTO emtr_all_info
      (job_id, method, rep, iter, root, ll_init, ll_final,
       root_prob_init, root_prob_final, trans_prob_init, trans_prob_final, ecd_ll_per_iter, raw_json)
    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
    ON DUPLICATE KEY UPDATE
      rep = VALUES(rep),
      iter = VALUES(iter),
      root = VALUES(root),
      ll_init = VALUES(ll_init),
      ll_final = VALUES(ll_final),
      root_prob_init = VALUES(root_prob_init),
      root_prob_final = VALUES(root_prob_final),
      trans_prob_init = VALUES(trans_prob_init),
      trans_prob_final = VALUES(trans_prob_final),
      ecd_ll_per_iter = VALUES(ecd_ll_per_iter),
      raw_json = VALUES(raw_json),
      updated_at = CURRENT_TIMESTAMP(6)
  `;

  await conn.execute(sql, [
    job_id, methodKey, rep, iter, root, ll_init, ll_final,
    JSON.stringify(root_prob_init),
    JSON.stringify(root_prob_final),
    JSON.stringify(trans_prob_init),
    JSON.stringify(trans_prob_final),
    JSON.stringify(ecd_ll_per_iter),
    JSON.stringify(raw_json),
  ]);
}

/* -------------------- POST: upsert best structs -------------------- */

export async function POST(req: NextRequest) {
  try {
    const ct = req.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      return NextResponse.json({ ok: false, error: "Content-Type must be application/json" }, { status: 415 });
    }

    const body: unknown = await req.json();
    const obj = asObj(body);
    if (!obj) return NextResponse.json({ ok: false, error: "Invalid JSON body" }, { status: 400 });

    const job_id = String(obj.job_id ?? "").trim();
    if (!job_id) return NextResponse.json({ ok: false, error: "job_id is required" }, { status: 400 });

    // Normalize any worker shape into {pars?, dirichlet?, hss?}
    const { pars, dirichlet, hss } = normalizeIncoming(obj);

    if (!pars && !dirichlet && !hss) {
      return NextResponse.json({ ok: false, error: "No best structs found (pars/dirichlet/hss missing)" }, { status: 400 });
    }

    const conn = db();
    if (pars)      await upsertOne(conn, job_id, "parsimony", { ...pars,      method: pars.method ?? "Parsimony" });
    if (dirichlet) await upsertOne(conn, job_id, "dirichlet", { ...dirichlet, method: dirichlet.method ?? "Dirichlet" });
    if (hss)       await upsertOne(conn, job_id, "hss",       { ...hss,       method: hss.method ?? "HSS" });

    return NextResponse.json({ ok: true }, { headers: { "cache-control": "no-store" } });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}

/* -------------------- GET: fetch best structs by job + Dirichlet params -------------------- */

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job_id = String(searchParams.get("job") ?? "").trim();
    if (!job_id) {
      return NextResponse.json({ ok: false, error: "missing ?job=" }, { status: 400 });
    }

    const conn = db();

    // 1) Dirichlet hyperparameters (independent numeric columns)
    const rsJob = (await conn.execute(
      `SELECT d_pi_1, d_pi_2, d_pi_3, d_pi_4,
              d_m_1,  d_m_2,  d_m_3,  d_m_4
         FROM emtr_jobs
        WHERE job_id = ?
        LIMIT 1`,
      [job_id]
    )) as ExecResult;

    const jobRows = rowsFromResult<JobParamRow>(rsJob, coerceJobParamRow);
    const jobParams = jobRows[0];
    const D_pi = jobParams
      ? [jobParams.d_pi_1, jobParams.d_pi_2, jobParams.d_pi_3, jobParams.d_pi_4].map((v) =>
          Number.isFinite(v as number) ? (v as number) : null
        )
      : [null, null, null, null];
    const D_M = jobParams
      ? [jobParams.d_m_1, jobParams.d_m_2, jobParams.d_m_3, jobParams.d_m_4].map((v) =>
          Number.isFinite(v as number) ? (v as number) : null
        )
      : [null, null, null, null];

    // 2) Fetch rows (all methods) for this job.
    // If multiple reps were written, you can change this to pick the best per method.
    const rsBest = (await conn.execute(
      `SELECT method, rep, iter, root, ll_init, ll_final,
              root_prob_init, root_prob_final,
              trans_prob_init, trans_prob_final,
              ecd_ll_per_iter, raw_json
         FROM emtr_all_info
        WHERE job_id = ?`,
      [job_id]
    )) as ExecResult;

    const rows = rowsFromResult<BestRow>(rsBest, coerceBestRow);

    let pars: EMStruct | null = null;
    let dirichlet: EMStruct | null = null;
    let hss: EMStruct | null = null;

    for (const row of rows) {
      const method = String(row.method ?? "").toLowerCase();

      // Prefer raw_json if present; else reconstruct from columns
      let em: EMStruct | null = null;
      if (row.raw_json) {
        try { em = JSON.parse(row.raw_json) as EMStruct; } catch { em = null; }
      }
      if (!em) {
        em = {
          method: row.method ?? undefined,
          rep: row.rep ?? null,
          iter: row.iter ?? null,
          ll_init: row.ll_init ?? null,
          ll_final: row.ll_final ?? null,
          root: row.root ?? null,
          root_prob_init: parseMaybeJson<number[]>(row.root_prob_init),
          root_prob_final: parseMaybeJson<number[]>(row.root_prob_final),
          trans_prob_init: parseMaybeJson<Record<string, unknown>>(row.trans_prob_init),
          trans_prob_final: parseMaybeJson<Record<string, unknown>>(row.trans_prob_final),
          ecd_ll_per_iter: parseMaybeJson<Record<string, number> | number[]>(row.ecd_ll_per_iter),
        };
      }

      if (method === "parsimony") pars = em;
      else if (method === "dirichlet") dirichlet = em;
      else if (method === "hss") hss = em;
    }

    return NextResponse.json(
      { ok: true, job_id, D_pi, D_M, pars, dirichlet, hss },
      { headers: { "cache-control": "no-store" } }
    );
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
