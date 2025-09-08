import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

// Run on Vercel Edge (compatible with PlanetScale serverless driver)
export const runtime = "edge";

type Body = {
  job_id: string;
  rep: number | null;
  iter?: number | null;
  root: string;

  aa_ll_init?: number | null;
  aa_ll_final?: number | null;

  aa_root_prob_init?: unknown;
  aa_root_prob_final?: unknown;

  aa_exchangeability_init?: unknown;
  aa_exchangeability_final?: unknown;

  aa_ll_per_iter?: unknown;
  raw_json?: unknown;
};

function asInt(x: unknown): number | null {
  const n = typeof x === "number" ? x : parseInt(String(x ?? ""), 10);
  return Number.isFinite(n) ? n : null;
}
function asNum(x: unknown): number | null {
  const n = typeof x === "number" ? x : Number(String(x ?? ""));
  return Number.isFinite(n) ? n : null;
}
function asStr(x: unknown): string {
  const s = String(x ?? "").trim();
  return s;
}
function asJsonText(x: unknown): string | null {
  if (x == null) return null;
  try { return JSON.stringify(x); } catch { return null; }
}

export async function POST(req: NextRequest) {
  try {
    const body = (await req.json()) as Body;

    const job_id = asStr(body.job_id);
    const rep = asInt(body.rep);
    const iter = asInt(body.iter);
    const root = asStr(body.root);

    const aa_ll_init  = asNum(body.aa_ll_init);
    const aa_ll_final = asNum(body.aa_ll_final);

  
    const aa_root_prob_init        = asJsonText(body.aa_root_prob_init);
    const aa_root_prob_final       = asJsonText(body.aa_root_prob_final);
    const aa_exchangeability_init  = asJsonText(body.aa_exchangeability_init);
    const aa_exchangeability_final = asJsonText(body.aa_exchangeability_final);
    const aa_ll_per_iter           = asJsonText(body.aa_ll_per_iter);
    const raw_json                 = asJsonText(body.raw_json ?? body);

  
    if (!job_id) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }
    if (!root) {
      return NextResponse.json({ ok: false, error: "Missing root" }, { status: 400 });
    }
    if (!Number.isFinite(rep as number)) {
      return NextResponse.json({ ok: false, error: "Missing or invalid rep" }, { status: 400 });
    }

    const conn = db();

    const sql = `
      INSERT INTO emtr_aa_info (
        job_id, rep, iter, root,
        aa_ll_init, aa_ll_final,
        aa_root_prob_init, aa_root_prob_final,
        aa_exchangeability_init, aa_exchangeability_final,
        aa_ll_per_iter, raw_json
      )
      VALUES (
        ?, ?, ?, ?,
        ?, ?,
        CAST(? AS JSON), CAST(? AS JSON),
        CAST(? AS JSON), CAST(? AS JSON),
        CAST(? AS JSON), CAST(? AS JSON)
      )
      ON DUPLICATE KEY UPDATE
        iter = VALUES(iter),
        aa_ll_init  = VALUES(aa_ll_init),
        aa_ll_final = VALUES(aa_ll_final),
        aa_root_prob_init        = VALUES(aa_root_prob_init),
        aa_root_prob_final       = VALUES(aa_root_prob_final),
        aa_exchangeability_init  = VALUES(aa_exchangeability_init),
        aa_exchangeability_final = VALUES(aa_exchangeability_final),
        aa_ll_per_iter = VALUES(aa_ll_per_iter),
        raw_json       = VALUES(raw_json),
        updated_at     = CURRENT_TIMESTAMP(6)
    `;

    const params: Array<string | number | null> = [
      job_id, rep, iter, root,
      aa_ll_init, aa_ll_final,
      aa_root_prob_init, aa_root_prob_final,
      aa_exchangeability_init, aa_exchangeability_final,
      aa_ll_per_iter, raw_json,
    ];

    await conn.execute(sql, params);

    return NextResponse.json({ ok: true });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : String(e);
    return NextResponse.json({ ok: false, error: msg.slice(0, 400) }, { status: 500 });
  }
}
