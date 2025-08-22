// app/api/best/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

/* -------------------- Types -------------------- */

type EMStruct = {
  method?: string;                                   // "Parsimony" | "Dirichlet" | "SSH"
  rep?: number | null;
  ecd_ll_per_iter?: Record<string, number> | number[] | null;
  ll_final?: number | null;
  root_name?: string | null;
  root_prob_init?: number[] | null;                  // length 4
  root_prob_final?: number[] | null;                 // length 4
  trans_prob_init?: Record<string, unknown> | null;  // keep as JSON
  trans_prob_final?: Record<string, unknown> | null; // keep as JSON
};

type BestRow = {
  method: string | null;
  rep: number | null;
  root_name: string | null;
  ll_final: number | null;
  root_prob_init: string | null;
  root_prob_final: string | null;
  trans_prob_init: string | null;
  trans_prob_final: string | null;
  ecd_ll_per_iter: string | null;
  raw_json: string | null;
};

/* -------------------- Helpers -------------------- */

function asObj(x: unknown): Record<string, unknown> | null {
  return x && typeof x === "object" ? (x as Record<string, unknown>) : null;
}

function coerceBestRow(u: unknown): BestRow {
  const o = (u && typeof u === "object" ? (u as Record<string, unknown>) : {}) as Record<
    string,
    unknown
  >;

  const toStr = (v: unknown): string | null =>
    v === null || v === undefined ? null : String(v);
  const toNum = (v: unknown): number | null =>
    v === null || v === undefined ? null : Number(v);

  return {
    method: toStr(o.method),
    rep: toNum(o.rep),
    root_name: toStr(o.root_name),
    ll_final: toNum(o.ll_final),
    root_prob_init: toStr(o.root_prob_init),
    root_prob_final: toStr(o.root_prob_final),
    trans_prob_init: toStr(o.trans_prob_init),
    trans_prob_final: toStr(o.trans_prob_final),
    ecd_ll_per_iter: toStr(o.ecd_ll_per_iter),
    raw_json: toStr(o.raw_json),
  };
}

function rowsFromResult<T>(
  rs: unknown,
  map: (u: unknown) => T
): T[] {
  const base = rs && typeof rs === "object" ? (rs as { rows?: unknown }) : null;
  if (!base || !Array.isArray(base.rows)) return [];
  return base.rows.map(map);
}

function parseMaybeJson<T>(v: unknown): T | null {
  if (v == null) return null;
  if (typeof v === "string") {
    try {
      return JSON.parse(v) as T;
    } catch {
      return null;
    }
  }
  // Already parsed (driver may return JSON columns as objects)
  return (v as T) ?? null;
}

/* -------------------- DB upsert -------------------- */

async function upsertOne(
  conn: ReturnType<typeof db>,
  job_id: string,
  methodKey: "parsimony" | "dirichlet" | "ssh",
  em: EMStruct
) {
  const rep = Number.isFinite(em.rep as number) ? (em.rep as number) : null;
  const root_name = typeof em.root_name === "string" ? em.root_name : null;
  const ll_final = Number.isFinite(em.ll_final as number) ? (em.ll_final as number) : null;

  const root_prob_init = em.root_prob_init ?? null;
  const root_prob_final = em.root_prob_final ?? null;
  const trans_prob_init = em.trans_prob_init ?? null;
  const trans_prob_final = em.trans_prob_final ?? null;
  const ecd_ll_per_iter = em.ecd_ll_per_iter ?? null;

  const raw_json = em; // store full struct as well

  await conn.execute(
    `INSERT INTO emtr_best
       (job_id, method, rep, root_name, ll_final,
        root_prob_init, root_prob_final, trans_prob_init, trans_prob_final, ecd_ll_per_iter, raw_json)
     VALUES (?,?,?,?,?,
             ?,?,?,?, ?,?)
     ON DUPLICATE KEY UPDATE
       rep=VALUES(rep),
       root_name=VALUES(root_name),
       ll_final=VALUES(ll_final),
       root_prob_init=VALUES(root_prob_init),
       root_prob_final=VALUES(root_prob_final),
       trans_prob_init=VALUES(trans_prob_init),
       trans_prob_final=VALUES(trans_prob_final),
       ecd_ll_per_iter=VALUES(ecd_ll_per_iter),
       raw_json=VALUES(raw_json),
       updated_at=CURRENT_TIMESTAMP(6)`,
    [
      job_id,
      methodKey,
      rep,
      root_name,
      ll_final,
      JSON.stringify(root_prob_init),
      JSON.stringify(root_prob_final),
      JSON.stringify(trans_prob_init),
      JSON.stringify(trans_prob_final),
      JSON.stringify(ecd_ll_per_iter),
      JSON.stringify(raw_json),
    ]
  );
}

/* -------------------- POST: upsert best structs -------------------- */

export async function POST(req: NextRequest) {
  try {
    const ct = req.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      return NextResponse.json(
        { ok: false, error: "Content-Type must be application/json" },
        { status: 415 }
      );
    }

    const body: unknown = await req.json();
    const obj = asObj(body);
    if (!obj) return NextResponse.json({ ok: false, error: "Invalid JSON body" }, { status: 400 });

    const job_id = String(obj.job_id ?? "").trim();
    if (!job_id) {
      return NextResponse.json({ ok: false, error: "job_id is required" }, { status: 400 });
    }

    const pars = asObj(obj.pars);
    const diri = asObj(obj.dirichlet);
    const ssh = asObj(obj.ssh);

    if (!pars && !diri && !ssh) {
      return NextResponse.json(
        { ok: false, error: "No best structs found (pars/dirichlet/ssh missing)" },
        { status: 400 }
      );
    }

    const conn = db();

    if (pars) {
      await upsertOne(conn, job_id, "parsimony", {
        ...pars,
        method: (pars.method as string) ?? "Parsimony",
      } as EMStruct);
    }
    if (diri) {
      await upsertOne(conn, job_id, "dirichlet", {
        ...diri,
        method: (diri.method as string) ?? "Dirichlet",
      } as EMStruct);
    }
    if (ssh) {
      await upsertOne(conn, job_id, "ssh", {
        ...ssh,
        method: (ssh.method as string) ?? "SSH",
      } as EMStruct);
    }

    return NextResponse.json({ ok: true }, { headers: { "cache-control": "no-store" } });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}

/* -------------------- GET: fetch best structs by job -------------------- */

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job_id = String(searchParams.get("job") ?? "").trim();
    if (!job_id) {
      return NextResponse.json({ ok: false, error: "missing ?job=" }, { status: 400 });
    }

    const conn = db();
    const rs = await conn.execute(
      `SELECT method, rep, root_name, ll_final,
              root_prob_init, root_prob_final,
              trans_prob_init, trans_prob_final,
              ecd_ll_per_iter, raw_json
         FROM emtr_best
        WHERE job_id = ?`,
      [job_id]
    );

    const rows = rowsFromResult<BestRow>(rs, coerceBestRow);

    let pars: EMStruct | null = null;
    let dirichlet: EMStruct | null = null;
    let ssh: EMStruct | null = null;

    for (const row of rows) {
      const method = String(row.method ?? "").toLowerCase();

      // Prefer the raw_json blob if present; otherwise reconstruct from JSON columns.
      let em: EMStruct | null = null;
      if (row.raw_json) {
        try {
          em = JSON.parse(row.raw_json) as EMStruct;
        } catch {
          em = null;
        }
      }
      if (!em) {
        em = {
          method: row.method ?? undefined,
          rep: row.rep ?? null,
          root_name: row.root_name ?? null,
          ll_final: row.ll_final ?? null,
          root_prob_init: parseMaybeJson<number[]>(row.root_prob_init),
          root_prob_final: parseMaybeJson<number[]>(row.root_prob_final),
          trans_prob_init: parseMaybeJson<Record<string, unknown>>(row.trans_prob_init),
          trans_prob_final: parseMaybeJson<Record<string, unknown>>(row.trans_prob_final),
          ecd_ll_per_iter: parseMaybeJson<Record<string, number> | number[]>(
            row.ecd_ll_per_iter
          ),
        };
      }

      if (method === "parsimony") pars = em;
      else if (method === "dirichlet") dirichlet = em;
      else if (method === "ssh") ssh = em;
    }

    return NextResponse.json(
      { ok: true, job_id, pars, dirichlet, ssh },
      { headers: { "cache-control": "no-store" } }
    );
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
