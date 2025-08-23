// app/api/ecdll/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Row = {
  ecd_ll_per_iter: string | unknown | null; // DB might store JSON text
};

type Pt = [number, number];

function normalizePoints(x: unknown): Pt[] {
  const src = typeof x === "string" ? safeParse(x) : x;
  if (!Array.isArray(src)) return [];
  const out: Pt[] = [];
  for (const row of src) {
    if (Array.isArray(row) && row.length >= 2) {
      const it = Number(row[0]);
      const ll = Number(row[1]);
      if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
    } else if (row && typeof row === "object") {
      const it = Number((row as Record<string, unknown>).iter);
      const ll = Number((row as Record<string, unknown>).ll);
      if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
    }
  }
  out.sort((a, b) => a[0] - b[0]);
  return out;
}

function safeParse(s: string): unknown {
  try { return JSON.parse(s); } catch { return s; }
}

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || "").trim();
    const method = (searchParams.get("method") || "dirichlet").trim();

    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }
    if (!method) {
      return NextResponse.json({ ok: false, error: "Missing method" }, { status: 400 });
    }

    const sql = `
      SELECT ecd_ll_per_iter
      FROM emtr_best_rep
      WHERE job_id = ? AND method = ?
      ORDER BY created_at DESC
      LIMIT 1
    `;

    const conn = db();
    const { rows } = await conn.execute<Row>(sql, [jobId, method]);

    if (!rows || rows.length === 0) {
      return NextResponse.json({ ok: true, job_id: jobId, method, points: [] });
    }

    const raw = rows[0]?.ecd_ll_per_iter ?? null;
    const points = normalizePoints(raw);

    return NextResponse.json({ ok: true, job_id: jobId, method, points });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
