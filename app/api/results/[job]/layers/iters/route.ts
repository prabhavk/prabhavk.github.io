// app/api/results/[job]/layers/iters/route.ts
import { NextResponse } from "next/server";
import { query as q } from "@/lib/db";

export const runtime = "nodejs";

type Row = {
  method: "parsimony" | "dirichlet" | "hss";
  layer: 0 | 1 | 2;
  iter: number | null;
};

export async function GET(req: Request) {
  try {
    const { pathname, search } = new URL(req.url);
    // /api/results/<job>/layers/iters
    const m = pathname.match(/\/api\/results\/([^/]+)\/layers\/iters\/?$/);
    const job = decodeURIComponent(m?.[1] ?? "");

    const u = new URL(req.url);
    const rep = Number(u.searchParams.get("rep"));
    const root = (u.searchParams.get("root") || "").trim();

    if (!job || !Number.isFinite(rep) || !root) {
      return NextResponse.json({ ok: false, error: "job, rep, root required" }, { status: 400 });
    }

    const rows = await q<Row>(
      `SELECT method, layer, iter
         FROM emtr_trifle_layers
        WHERE job_id = ? AND rep = ? AND root = ?
        ORDER BY method, layer`,
      [job, rep, root]
    );

    return NextResponse.json({ ok: true, rows: rows ?? [] });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
