// app/api/pattern-weights/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type Body = {
  job_id?: string;
  num_patterns?: number;
  weights?: number[];
  cum_weights?: number[];
};

function errMsg(e: unknown): string {
  if (e instanceof Error) return e.message;
  if (typeof e === "string") return e;
  try { return JSON.stringify(e); } catch { return "Unknown error"; }
}
function pickStr(...vals: Array<unknown>): string {
  for (const v of vals) if (typeof v === "string" && v.trim()) return v.trim();
  return "";
}
function isNumArray(x: unknown): x is number[] {
  return Array.isArray(x) && x.every(v => Number.isFinite(Number(v)));
}

// ---------- POST: upsert full arrays per job ----------
export async function POST(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobFromQuery  = (searchParams.get("job_id") ?? searchParams.get("job") ?? "").trim();
    const jobFromHeader = (req.headers.get("x-job-id") ?? "").trim();

    const body = (await req.json().catch(() => null)) as Body | null;
    if (!body) return NextResponse.json({ ok: false, error: "Invalid JSON" }, { status: 400 });

    const job_id = pickStr(jobFromQuery, jobFromHeader, body.job_id);
    if (!job_id) return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });

    const weights = body.weights;
    const cums    = body.cum_weights;
    if (!isNumArray(weights) || !isNumArray(cums) || weights.length !== cums.length) {
      return NextResponse.json({ ok: false, error: "Expected numeric arrays: weights[] and cum_weights[] of equal length" }, { status: 400 });
    }

    const num_patterns = Number.isFinite(body.num_patterns) ? Number(body.num_patterns) : weights.length;

    // Store as JSON strings (works whether columns are JSON or TEXT)
    const weightsJson = JSON.stringify(weights.map(Number));
    const cumsJson    = JSON.stringify(cums.map(Number));

    const sql = `
      INSERT INTO emtr_pattern_weights (job_id, num_patterns, weights, cum_weights)
      VALUES (?, ?, ?, ?)
      ON DUPLICATE KEY UPDATE
        num_patterns = VALUES(num_patterns),
        weights      = VALUES(weights),
        cum_weights  = VALUES(cum_weights)
    `;
    await query(sql, [job_id, num_patterns, weightsJson, cumsJson]);

    return NextResponse.json({ ok: true, job_id, num_patterns, inserted: 1 });
  } catch (e) {
    return NextResponse.json({ ok: false, error: errMsg(e) }, { status: 400 });
  }
}

// ---------- GET: latest record for a job ----------
export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job_id = (searchParams.get("job_id") ?? searchParams.get("job") ?? "").trim();
    if (!job_id) return NextResponse.json({ error: "Missing job_id" }, { status: 400 });

    const rows = await query<{
      id: number; job_id: string; num_patterns: number;
      weights: string | unknown; cum_weights: string | unknown;
      created_at?: string; updated_at?: string;
    }>(
      `
      SELECT id, job_id, num_patterns, weights, cum_weights, created_at, updated_at
        FROM emtr_pattern_weights
       WHERE job_id = ?
       ORDER BY updated_at DESC, id DESC
       LIMIT 1
      `,
      [job_id]
    );

    if (!rows.length) return NextResponse.json({ error: "Not found" }, { status: 404 });

    const r = rows[0];
    const weights = Array.isArray(r.weights) ? r.weights : JSON.parse(String(r.weights || "[]"));
    const cum     = Array.isArray(r.cum_weights) ? r.cum_weights : JSON.parse(String(r.cum_weights || "[]"));

    return NextResponse.json({
      job_id: r.job_id,
      num_patterns: r.num_patterns,
      weights,
      cum_weights: cum,
      id: r.id,
      created_at: r.created_at,
      updated_at: r.updated_at,
    });
  } catch (e) {
    return NextResponse.json({ error: errMsg(e) }, { status: 500 });
  }
}
