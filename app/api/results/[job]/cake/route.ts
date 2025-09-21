// app/api/results/[job]/cake/route.ts
import { NextResponse } from "next/server";
import { query as q } from "@/lib/db";

export const runtime = "nodejs";

type M = "parsimony" | "dirichlet" | "hss";
type Row = { method: M; layer: number; ll_final: number | string | null };

const LAYER_KEY: Record<number, "coarse" | "medium" | "fine"> = {
  0: "coarse", 1: "medium", 2: "fine",
};

export async function GET(req: Request) {
  // Parse job from the pathname: /api/results/<job>/cake
  const { pathname, search } = new URL(req.url);
  const m = pathname.match(/\/api\/results\/([^/]+)\/cake\/?$/);
  const job = decodeURIComponent(m?.[1] ?? "");

  const u = new URL(req.url);
  const rep = Number(u.searchParams.get("rep"));
  const root = (u.searchParams.get("root") || "").trim();

  if (!job || !Number.isFinite(rep) || !root) {
    return NextResponse.json({
      ok: false,
      job_id: job || null,
      rep: Number.isFinite(rep) ? rep : null,
      root: root || null,
      pattern_counts: { coarse: 30, medium: 60, fine: 100 },
      ll: {},
    });
  }

  const rows = await q<Row>(
    `SELECT method, layer, MAX(ll_final) AS ll_final
       FROM emtr_trifle_layers
      WHERE job_id = ? AND rep = ? AND root = ?
      GROUP BY method, layer`,
    [job, rep, root]
  );

  const init = () => ({ coarse: null as number | null, medium: null as number | null, fine: null as number | null });
  const ll: Record<M, ReturnType<typeof init>> = {
    parsimony: init(), dirichlet: init(), hss: init(),
  };

  for (const r of rows ?? []) {
    const key = LAYER_KEY[r.layer];
    if (!key) continue;
    const v = r.ll_final == null ? null : Number(r.ll_final);
    if (Number.isFinite(v as number)) ll[r.method][key] = v as number;
  }

  // TODO: replace with real counts if/when you store them
  const pattern_counts = { coarse: 30, medium: 60, fine: 100 };

  return NextResponse.json({ ok: true, job_id: job, rep, root, pattern_counts, ll });
}
