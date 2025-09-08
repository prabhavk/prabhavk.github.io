// app/api/ecdll/roots/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Row = { root: string | null; cnt?: string | number | null };

// Normalize root strings; keep only non-empty strings.
function toRoot(x: unknown): string | null {
  if (typeof x !== "string") return null;
  const s = x.trim();
  return s.length ? s : null;
}

// Natural sort: h_# (internal nodes) first by numeric suffix, then others A→Z.
function naturalRootCompare(a: string, b: string): number {
  const ha = /^h_(\d+)$/.exec(a);
  const hb = /^h_(\d+)$/.exec(b);
  if (ha && hb) {
    const na = Number(ha[1]);
    const nb = Number(hb[1]);
    return na - nb;
  }
  if (ha) return -1; // a is internal, b not → a first
  if (hb) return 1;  // b is internal, a not → b first
  // both non-internal → case-insensitive lexical
  return a.toLowerCase().localeCompare(b.toLowerCase()) || a.localeCompare(b);
}

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const jobId =
      (searchParams.get("job_id") || searchParams.get("job") || "").trim();
    const method = (searchParams.get("method") || "").trim().toLowerCase(); // optional

    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }

    const conn = db();

    // Build SQL; include method filter only if provided. We also fetch counts so we can
    // pick a reasonable default (most frequent root) if you want to use it.
    let sql = `
      SELECT root, COUNT(*) AS cnt
      FROM emtr_all_info
      WHERE job_id = ?
        AND root IS NOT NULL
        AND root <> ''
    `;
    const params: Array<string | number> = [jobId];

    if (method) {
      sql += ` AND LOWER(method) = ?`;
      params.push(method);
    }

    sql += ` GROUP BY root`;

    const { rows } = await conn.execute<Row>(sql, params);

    const roots = (rows ?? [])
      .map((r) => toRoot(r.root))
      .filter((v): v is string => v !== null);

    // Dedup (defensive) and natural sort (h_# first)
    const uniq = Array.from(new Set(roots));
    uniq.sort(naturalRootCompare);

    // Choose a preferred default root:
    //   1) If we have counts, prefer the most frequent root,
    //   2) else fall back to the first after sorting.
    let preferred: string | null = null;
    if (rows && rows.length) {
      // Find row with max count among valid roots
      let best: { root: string; cnt: number } | null = null;
      for (const r of rows) {
        const name = toRoot(r.root);
        if (!name) continue;
        const c = Number(r.cnt);
        const cnt = Number.isFinite(c) ? c : 0;
        if (!best || cnt > best.cnt) best = { root: name, cnt };
      }
      preferred = best?.root ?? null;
    }
    if (!preferred) preferred = uniq[0] ?? null;

    return NextResponse.json({
      ok: true,
      job_id: jobId,
      method: method || null,
      roots: uniq,
      root: preferred,
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
