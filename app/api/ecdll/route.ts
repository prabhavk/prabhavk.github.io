// app/api/ecdll/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Row = {
  ecd_ll_per_iter: string | null;
  root_name: string | null;
  root_prob_final: string | null;
  ll_init: number | null;
  ll_final: number | null;
  rep: number | null;
};

function parsePoints(jsonish: string | null): [number, number][] {
  if (!jsonish) return [];
  try {
    const raw = JSON.parse(jsonish);
    if (!Array.isArray(raw)) return [];
    const out: [number, number][] = [];
    for (const r of raw) {
      if (Array.isArray(r) && r.length >= 2) {
        const it = Number(r[0]);
        const ll = Number(r[1]);
        if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
      } else if (r && typeof r === "object") {
        const obj = r as Record<string, unknown>;
        const it = Number(obj["iter"] as number | string);
        const ll = Number(obj["ll"] as number | string);
        if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
      }
    }
    out.sort((a, b) => a[0] - b[0]);
    return out;
  } catch {
    return [];
  }
}



export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || "").trim();
    const method = (searchParams.get("method") || "").trim().toLowerCase();
    const repStr = searchParams.get("rep");
    const rep = repStr ? Number(repStr) : null;

    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }
    if (!method) {
      return NextResponse.json({ ok: false, error: "Missing method" }, { status: 400 });
    }

    // Pull a single run (row) for this job/method/rep.
    // We assume you now store **one row per run** in emtr_best_rep with a rep column.
    const sql = `
      SELECT
        ecd_ll_per_iter,
        root_name,
        root_prob_final,
        ll_init,
        ll_final,
        rep
      FROM emtr_best_rep
      WHERE job_id = ? AND method = ? ${rep !== null && Number.isFinite(rep) ? "AND rep = ?" : ""}
      ORDER BY created_at DESC
      LIMIT 1
    `;

    const conn = db();
    const { rows } = await conn.execute<Row>(
      sql,
      rep !== null && Number.isFinite(rep) ? [jobId, method, rep] : [jobId, method]
    );

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
      rep: row.rep,
      points,
      summary: {
        root: row.root_name ?? null,
        num_iterations: points.length,
        ll_init: row.ll_init,
        ecd_ll_first: ecd_first,
        ecd_ll_final: ecd_final,
        ll_final: row.ll_final,
        root_prob_final: row.root_prob_final ? JSON.parse(row.root_prob_final) : null,
      },
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
