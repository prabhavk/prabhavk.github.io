// app/api/trifle/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

/* ------------------------------- Types ------------------------------- */

type Method = "parsimony" | "dirichlet" | "hss";

type Layer = {
  layer: 0 | 1 | 2;
  iter: number | null;
  ll_initial: number | null;
  ll_final: number | null;
  root_prob_final?: unknown | null;
  trans_prob_final?: unknown | null;
  ecd_ll_per_iter?: Record<string, number> | null;
};

/* --------------------------- Type utilities -------------------------- */

function isRecord(v: unknown): v is Record<string, unknown> {
  return typeof v === "object" && v !== null;
}

function safeParseJSON(s: string): unknown | null {
  try { return JSON.parse(s); } catch { return null; }
}

function toInt(x: unknown): number | null {
  if (x === null || x === undefined || x === "") return null;
  const n = Number.parseInt(String(x), 10);
  return Number.isFinite(n) ? n : null;
}

function toNum(x: unknown): number | null {
  if (x === null || x === undefined || x === "") return null;
  const n = Number(x);
  return Number.isFinite(n) ? n : null;
}

function isMethod(v: unknown): v is Method {
  return v === "parsimony" || v === "dirichlet" || v === "hss";
}

/* -------------------------- Layer normalization ---------------------- */

function coerceLayers(body: Record<string, unknown>): Layer[] | null {
  // Case A: preferred layers[]
  if (Array.isArray(body["layers"])) {
    const out: Layer[] = [];
    for (const raw of body["layers"] as unknown[]) {
      if (!isRecord(raw)) continue;
      const layer = toInt(raw["layer"]);
      if (layer !== 0 && layer !== 1 && layer !== 2) continue;
      out.push({
        layer,
        iter: toInt(raw["iter"]),
        ll_initial: toNum(raw["ll_initial"]),
        ll_final: toNum(raw["ll_final"]),
        root_prob_final: raw["root_prob_final"] ?? null,
        trans_prob_final: raw["trans_prob_final"] ?? null,
        ecd_ll_per_iter: isRecord(raw["ecd_ll_per_iter"])
          ? (raw["ecd_ll_per_iter"] as Record<string, number>)
          : null,
      });
    }
    return out.length ? out : null;
  }

  // Case B: legacy vector-per-field
  const iterV = Array.isArray(body["iter_trifle"]) ? (body["iter_trifle"] as unknown[]) : [];
  const lliV  = Array.isArray(body["ll_initial_trifle"]) ? (body["ll_initial_trifle"] as unknown[]) : [];
  const llfV  = Array.isArray(body["ll_final_trifle"]) ? (body["ll_final_trifle"] as unknown[]) : [];
  const rpfV  = Array.isArray(body["root_prob_final_trifle"]) ? (body["root_prob_final_trifle"] as unknown[]) : [];
  const tpfV  = Array.isArray(body["trans_prob_final_trifle"]) ? (body["trans_prob_final_trifle"] as unknown[]) : [];
  const ecdV  = Array.isArray(body["ecd_ll_per_iter_for_trifle"]) ? (body["ecd_ll_per_iter_for_trifle"] as unknown[]) : [];

  const out: Layer[] = [];
  for (let k = 0 as 0 | 1 | 2; k < 3; k = (k + 1) as 0 | 1 | 2) {
    const iter = toInt(iterV[k]);
    const lli  = toNum(lliV[k]);
    const llf  = toNum(llfV[k]);

    let ecd: Record<string, number> | null = null;
    const raw = ecdV[k];
    if (isRecord(raw)) {
      ecd = raw as Record<string, number>;
    } else if (Array.isArray(raw)) {
      const obj: Record<string, number> = {};
      for (const p of raw) {
        if (Array.isArray(p) && p.length >= 2) {
          const ki = toInt(p[0]); const vv = toNum(p[1]);
          if (ki !== null && vv !== null) obj[String(ki)] = vv;
        }
      }
      ecd = Object.keys(obj).length ? obj : null;
    }

    if (iter !== null || lli !== null || llf !== null || rpfV[k] != null || tpfV[k] != null || ecd != null) {
      out.push({
        layer: k,
        iter,
        ll_initial: lli,
        ll_final: llf,
        root_prob_final: (rpfV[k] ?? null) as unknown,
        trans_prob_final: (tpfV[k] ?? null) as unknown,
        ecd_ll_per_iter: ecd,
      });
    }
  }
  return out.length ? out : null;
}

/* --------------------------------- POST -------------------------------- */

export async function POST(req: NextRequest) {
  try {
    // Parse body safely
    const bodyText = await req.text().catch(() => "");
    const parsed = bodyText ? safeParseJSON(bodyText) : {};
    const b: Record<string, unknown> = isRecord(parsed) ? parsed : {};

    // job_id: prefer header, then body
    const jobId =
      req.headers.get("x-job-id") ??
      (typeof b["job_id"] === "string" ? (b["job_id"] as string) : null);

    if (!jobId || !jobId.trim()) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }

    const methodRaw = b["method"];
    if (!isMethod(methodRaw)) {
      return NextResponse.json({ ok: false, error: "Invalid method" }, { status: 400 });
    }
    const method: Method = methodRaw;

    const rep = toInt(b["rep"]);
    const root = typeof b["root"] === "string" ? (b["root"] as string).trim() : "";

    if (rep === null || !root) {
      return NextResponse.json({ ok: false, error: "Missing rep/root" }, { status: 400 });
    }

    // Optional per-rep initial fields
    const ll_init = toNum(b["ll_init"]);
    const root_prob_init = b["root_prob_init"] ?? null;
    const trans_prob_init = b["trans_prob_init"] ?? null;

    // Layers
    const layers = coerceLayers(b);
    if (!layers || layers.length === 0) {
      return NextResponse.json({ ok: false, error: "Missing layers" }, { status: 400 });
    }

    // ---------- Upsert header ----------
    await query(
      `INSERT INTO emtr_trifle_runs
         (job_id, method, rep, root, ll_init, root_prob_init, trans_prob_init)
       VALUES
         (?, ?, ?, ?, ?, CAST(? AS JSON), CAST(? AS JSON))
       ON DUPLICATE KEY UPDATE
         ll_init         = VALUES(ll_init),
         root_prob_init  = VALUES(root_prob_init),
         trans_prob_init = VALUES(trans_prob_init)`,
      [
        jobId,
        method,
        rep,
        root,
        ll_init,
        JSON.stringify(root_prob_init),
        JSON.stringify(trans_prob_init),
      ]
    );

    // ---------- Upsert layers ----------
    let layersUpserted = 0;
    for (const L of layers) {
      try {
        await query(
          `INSERT INTO emtr_trifle_layers
             (job_id, method, rep, root, layer,
              iter, ll_initial, ll_final, root_prob_final, trans_prob_final, ecd_ll_per_iter)
           VALUES
             (?, ?, ?, ?, ?, ?, ?, ?, CAST(? AS JSON), CAST(? AS JSON), CAST(? AS JSON))
           ON DUPLICATE KEY UPDATE
             iter             = VALUES(iter),
             ll_initial       = VALUES(ll_initial),
             ll_final         = VALUES(ll_final),
             root_prob_final  = VALUES(root_prob_final),
             trans_prob_final = VALUES(trans_prob_final),
             ecd_ll_per_iter  = VALUES(ecd_ll_per_iter)`,
          [
            jobId,
            method,
            rep,
            root,
            L.layer,
            L.iter,
            L.ll_initial,
            L.ll_final,
            JSON.stringify(L.root_prob_final ?? null),
            JSON.stringify(L.trans_prob_final ?? null),
            JSON.stringify(L.ecd_ll_per_iter ?? null),
          ]
        );
        layersUpserted++;
      } catch (e) {
        // If the layers table/columns aren't present yet, skip gracefully.
        const msg = String((e as Error)?.message ?? e);
        if (
          msg.includes("doesn't exist") ||
          msg.includes("Unknown table") ||
          msg.includes("Unknown column") ||
          msg.includes("no such table")
        ) {
          break;
        }
        throw e;
      }
    }

    return NextResponse.json({
      ok: true,
      data: {
        job_id: jobId,
        method,
        rep,
        root,
        ll_init,
        layers_received: layers.length,
        layers_upserted: layersUpserted,
      },
    });
  } catch (e) {
    const msg = String((e as Error)?.message ?? e);
    return NextResponse.json({ ok: false, error: `DB error: ${msg}` }, { status: 500 });
  }
}
