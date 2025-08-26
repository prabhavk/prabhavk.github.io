// app/api/ecdll/reps/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

// PlanetScale/MySQL often returns numbers as strings; handle both.
type Row = { rep: string | number | null };

// Safe integer coercion (returns null if not finite)
function toInt(x: unknown): number | null {
  const n = Number(x);
  return Number.isFinite(n) ? Math.trunc(n) : null;
}

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId = (searchParams.get("job_id") || searchParams.get("job") || "").trim();
    const method = (searchParams.get("method") || "").trim().toLowerCase(); // optional

    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }

    const conn = db();

    // Build SQL; include method filter only if provided
    let sql = `
      SELECT DISTINCT rep
      FROM emtr_all_info
      WHERE job_id = ?
        AND rep IS NOT NULL
    `;
    const params: Array<string | number> = [jobId];

    if (method) {
      sql += ` AND LOWER(method) = ?`;
      params.push(method);
    }

    sql += ` ORDER BY rep ASC`;

    const { rows } = await conn.execute<Row>(sql, params);

    // Coerce to integers; keep only finite ints
    const reps = (rows ?? [])
      .map((r) => toInt(r.rep))
      .filter((v): v is number => v !== null);

    // rows are DISTINCT + ORDER BY already, but normalize just in case
    const uniqSorted = Array.from(new Set(reps)).sort((a, b) => a - b);

    return NextResponse.json({
      ok: true,
      job_id: jobId,
      method: method || null,
      reps: uniqSorted,
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
