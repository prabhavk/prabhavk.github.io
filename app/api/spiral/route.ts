// app/api/spiral/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type MethodName = "Parsimony" | "Dirichlet" | "SSH";

// What the page expects
type FlowerRow = {
  root: string;
  rep: number;
  ll_init: number | null;
  ecd_ll_first: number | null;
  ecd_ll_final: number | null;
  ll_final: number | null;
};

type ApiFlower = {
  job_id: string;
  method: MethodName;
  rows: FlowerRow[];
  reps: number[];
};

function normalizeMethod(s: string): MethodName | null {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "Parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "Dirichlet";
  if (m.includes("ssh")) return "SSH";
  return null;
}

export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job = searchParams.get("job") ?? "";
    const methodInput = searchParams.get("method") ?? "";

    if (!job) return NextResponse.json({ error: "Missing ?job=<job_id>" }, { status: 400 });

    const norm = normalizeMethod(methodInput);
    if (!norm) return NextResponse.json({ error: "Invalid ?method" }, { status: 400 });

    // Match DB 'method' case-insensitively
    const methodLike = norm.toLowerCase(); // parsimony | dirichlet | ssh

    // Pull all reps/rows for job+method; only from completed jobs
    const rows = await query<FlowerRow>(
      `
      SELECT r.root,
             r.rep,
             r.ll_init,
             r.ecd_ll_first,
             r.ecd_ll_final,
             r.ll_final
        FROM emtr_rows r
        JOIN emtr_jobs j ON j.job_id = r.job_id
       WHERE r.job_id = ?
         AND LOWER(r.method) = ?
         AND j.status = 'completed'
       ORDER BY r.root, r.rep
      `,
      [job, methodLike]
    );

    // Distinct reps present
    const repSet = new Set<number>();
    for (const r of rows) { if (Number.isFinite(r.rep)) repSet.add(r.rep); }
    const reps = Array.from(repSet).sort((a, b) => a - b);

    const payload: ApiFlower = {
      job_id: job,
      method: norm,
      rows,
      reps,
    };

    return NextResponse.json(payload);
  } catch (err: unknown) {
    const message = err instanceof Error ? err.message : "Internal error";
    return NextResponse.json({ error: message }, { status: 500 });
  }
}
