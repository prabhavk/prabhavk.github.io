// app/api/trifle/layers/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { q } from "@/lib/db";

export const runtime = "nodejs";

// GET /api/trifle/layers?job_id=J123[&rep=1][&root=A]
export async function GET(req: NextRequest) {
  const url = new URL(req.url);
  const jobId = (url.searchParams.get("job_id") || "").trim();
  const repStr = url.searchParams.get("rep");
  const root = url.searchParams.get("root") || undefined;

  if (!jobId) {
    return NextResponse.json({ ok: false, error: "job_id is required" }, { status: 400 });
  }

  const args: (string | number)[] = [jobId];
  const where: string[] = ["job_id = ?"];
  if (repStr != null && repStr !== "") {
    const rep = Number(repStr);
    if (!Number.isInteger(rep)) {
      return NextResponse.json({ ok: false, error: "rep must be an integer" }, { status: 400 });
    }
    where.push("rep = ?");
    args.push(rep);
  }
  if (root) {
    where.push("root = ?");
    args.push(root);
  }

  type Row = {
    method: "parsimony" | "dirichlet" | "hss";
    rep: number;
    root: string;
    layer: 0 | 1 | 2;
    iter: number | null;
    ll_final: number | null;
  };

  const rows = await q<Row>(
    `SELECT method, rep, root, layer, iter, ll_final
       FROM emtr_trifle_layers
      WHERE ${where.join(" AND ")}
      ORDER BY method, rep, root, layer`,
    args
  );

  // Distinct reps/roots for selectors
  const reps = Array.from(new Set(rows.map(r => r.rep))).sort((a, b) => a - b);
  const roots = Array.from(new Set(rows.map(r => r.root))).sort();

  // Pull pattern weights (coarse, medium, fine) from emtr_jobs.pattern_counts JSON, if present
  // Falls back to {26, 42, 100} if not available.
  type PatRow = { coarse: number | null; medium: number | null; fine: number | null };
  let pat: PatRow[] = [];
  try {
    pat = await q<PatRow>(
      `SELECT
         JSON_EXTRACT(pattern_counts, '$.coarse') AS coarse,
         JSON_EXTRACT(pattern_counts, '$.medium') AS medium,
         JSON_EXTRACT(pattern_counts, '$.fine')   AS fine
       FROM emtr_jobs
       WHERE job_id = ?
       LIMIT 1`,
      [jobId]
    );
  } catch {
    // Column might not exist in your schema; ignore.
  }
  const pct = pat[0] || { coarse: null, medium: null, fine: null };
  const patternCounts = {
    coarse: Number.isFinite(Number(pct.coarse)) ? Number(pct.coarse) : 26,
    medium: Number.isFinite(Number(pct.medium)) ? Number(pct.medium) : 42,
    fine:   Number.isFinite(Number(pct.fine))   ? Number(pct.fine)   : 100,
  };

  return NextResponse.json({
    ok: true,
    job_id: jobId,
    patternCounts,
    reps,
    roots,
    rows,
  });
}
