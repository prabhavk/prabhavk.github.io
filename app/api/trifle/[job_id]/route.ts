import { NextResponse } from "next/server";
import mysql, { type ResultSetHeader, type SslOptions } from "mysql2/promise";

export const runtime = "nodejs";

/* ---------------------------- Types & guards ---------------------------- */

type Method = "parsimony" | "dirichlet" | "hss";

interface LayerIn {
  layer: 0 | 1 | 2;
  iter?: number | null;
  ll_initial?: number | null;   // ← added
  ll_final?: number | null;
  root_prob_final?: unknown;
  trans_prob_final?: unknown;
  ecd_ll_per_iter?: unknown;    // object map or array-of-pairs is OK; stored as JSON
}

interface TrifleBody {
  method: Method;
  root: string;
  rep: number;
  ll_init?: number | null;
  root_prob_init?: unknown;
  trans_prob_init?: unknown;
  layers: LayerIn[];
}

function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isMethod(v: unknown): v is Method {
  return v === "parsimony" || v === "dirichlet" || v === "hss";
}
function isLayer(x: unknown): x is LayerIn {
  if (!isObject(x)) return false;
  const layer = Number((x as { layer?: unknown }).layer);
  return layer === 0 || layer === 1 || layer === 2;
}
function isLayerArray(x: unknown): x is LayerIn[] {
  return Array.isArray(x) && x.every(isLayer);
}
const numOrNull = (v: unknown): number | null => {
  const n = Number(v);
  return Number.isFinite(n) ? n : null;
};

/* ---- robust param extraction (works for edge/node and promise params) --- */
function isPromise<T = unknown>(x: unknown): x is Promise<T> {
  return isObject(x) && "then" in x && typeof (x as { then?: unknown }).then === "function";
}
async function getJobIdFromContext(ctx: unknown): Promise<string | null> {
  if (!isObject(ctx) || !("params" in ctx)) return null;
  const p = (ctx as { params: unknown }).params;
  const params = isPromise(p) ? await p : p;
  if (!isObject(params)) return null;
  const jid = (params as Record<string, unknown>).job_id;
  return typeof jid === "string" ? jid : null;
}

/* ----------------------------- DB connection ---------------------------- */

function parseMySQLURL(urlStr: string): {
  host: string; port: number; user: string; password: string; database: string; ssl?: SslOptions;
} {
  const u = new URL(urlStr);
  const host = u.hostname;
  const port = u.port ? Number(u.port) : 3306;
  const user = decodeURIComponent(u.username);
  const password = decodeURIComponent(u.password);
  const database = u.pathname.replace(/^\//, "");
  let ssl: SslOptions | undefined;
  const sslParam = u.searchParams.get("ssl");
  if (sslParam) {
    try {
      const parsed = JSON.parse(sslParam) as Record<string, unknown>;
      const ra = typeof parsed.rejectUnauthorized === "boolean" ? parsed.rejectUnauthorized : true;
      ssl = { rejectUnauthorized: ra };
    } catch {}
  }
  return { host, port, user, password, database, ssl };
}

const DATABASE_URL = process.env.DATABASE_URL;
if (!DATABASE_URL) throw new Error("Missing DATABASE_URL");
const cfg = parseMySQLURL(DATABASE_URL);

const pool = mysql.createPool({
  host: cfg.host, port: cfg.port, user: cfg.user, password: cfg.password, database: cfg.database,
  ssl: cfg.ssl ?? { rejectUnauthorized: true }, waitForConnections: true, connectionLimit: 5,
});

/* --------------------------------- POST --------------------------------- */

export async function POST(req: Request, context: unknown) {
  const url = new URL(req.url);

  // prefer the `[job_id]` path param, then header, then `?job_id`
  const jobFromParam  = (await getJobIdFromContext(context)) ?? "";
  const jobFromHeader = (req.headers.get("x-job-id") ?? "").trim();
  const jobFromQuery  = (url.searchParams.get("job_id") ?? "").trim();
  const job_id = jobFromParam.trim() || jobFromHeader || jobFromQuery;

  // parse body safely
  const raw: unknown = await req.json().catch(() => null);
  if (!raw || !isObject(raw)) {
    return NextResponse.json({ ok: false, error: "Invalid JSON" }, { status: 400 });
  }

  // shape + validate
  const methodStr = String((raw as Record<string, unknown>).method ?? "").trim();
  const root      = String((raw as Record<string, unknown>).root ?? "").trim();
  const repNum    = Number((raw as Record<string, unknown>).rep);
  const layers    = isLayerArray((raw as Record<string, unknown>).layers)
    ? ((raw as Record<string, unknown>).layers as LayerIn[])
    : [];

  if (!job_id || !isMethod(methodStr) || !root || !Number.isFinite(repNum) || layers.length === 0) {
    return NextResponse.json({ ok: false, error: "Missing required fields" }, { status: 400 });
  }

  const ll_init         = numOrNull((raw as Record<string, unknown>).ll_init);
  const root_prob_init  = (raw as Record<string, unknown>).root_prob_init ?? null;
  const trans_prob_init = (raw as Record<string, unknown>).trans_prob_init ?? null;

  const conn = await pool.getConnection();
  try {
    await conn.beginTransaction();

    // Upsert run header and capture run_id (LAST_INSERT_ID trick)
    const [r] = await conn.query<ResultSetHeader>(
      `
      INSERT INTO emtr_trifle_runs
        (job_id, method, root, rep, ll_init, root_prob_init, trans_prob_init)
      VALUES (?,?,?,?,?,CAST(? AS JSON),CAST(? AS JSON))
      ON DUPLICATE KEY UPDATE
        ll_init = VALUES(ll_init),
        root_prob_init = VALUES(root_prob_init),
        trans_prob_init = VALUES(trans_prob_init),
        id = LAST_INSERT_ID(id)
      `,
      [
        job_id,
        methodStr,
        root,
        repNum,
        ll_init,
        root_prob_init ? JSON.stringify(root_prob_init) : null,
        trans_prob_init ? JSON.stringify(trans_prob_init) : null,
      ]
    );
    const runId = r.insertId;

    // Upsert each layer (now includes ll_initial)
    const insertLayerSql = `
      INSERT INTO emtr_trifle_layers
        (run_id, job_id, method, root, rep, layer,
         iter, ll_initial, ll_final, root_prob_final, trans_prob_final, ecd_ll_per_iter)
      VALUES (?,?,?,?,?,?,?,?,?,CAST(? AS JSON),CAST(? AS JSON),CAST(? AS JSON))
      ON DUPLICATE KEY UPDATE
        iter = VALUES(iter),
        ll_initial = VALUES(ll_initial),
        ll_final = VALUES(ll_final),
        root_prob_final = VALUES(root_prob_final),
        trans_prob_final = VALUES(trans_prob_final),
        ecd_ll_per_iter = VALUES(ecd_ll_per_iter),
        updated_at = CURRENT_TIMESTAMP(6)
    `;

    for (const L of layers) {
      const layerIdx = Number(L.layer);
      if (layerIdx !== 0 && layerIdx !== 1 && layerIdx !== 2) continue;

      await conn.query<ResultSetHeader>(insertLayerSql, [
        runId,
        job_id,
        methodStr,
        root,
        repNum,
        layerIdx,
        numOrNull(L.iter),
        numOrNull(L.ll_initial),  // ← write ll_initial
        numOrNull(L.ll_final),
        L.root_prob_final ? JSON.stringify(L.root_prob_final) : null,
        L.trans_prob_final ? JSON.stringify(L.trans_prob_final) : null,
        L.ecd_ll_per_iter ? JSON.stringify(L.ecd_ll_per_iter) : null,
      ]);
    }

    await conn.commit();
    return NextResponse.json({ ok: true, run_id: runId });
  } catch (e: unknown) {
    try { await conn.rollback(); } catch {}
    const msg = e instanceof Error ? e.message : "db error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  } finally {
    conn.release();
  }
}
