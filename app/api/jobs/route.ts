// app/api/jobs/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Num = number;

interface JobPayload {
  job_id: string;
  thr: Num;
  reps: Num;
  max_iter: Num;
  D_pi?: Num[]; // length 4 if present
  D_M?: Num[];  // length 4 if present
}

function asNumberArray(x: unknown): number[] | null {
  if (!Array.isArray(x)) return null;
  const arr = x.map((v) => Number(v));
  return arr.every((v) => Number.isFinite(v)) ? arr : null;
}

export async function POST(req: NextRequest) {
  try {
    const ct = req.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      return NextResponse.json(
        { ok: false, error: "Content-Type must be application/json" },
        { status: 415 }
      );
    }

    const raw = (await req.json()) as Partial<JobPayload>;
    const job_id = String(raw.job_id ?? "").trim();
    const thr = Number(raw.thr);
    const reps = Number(raw.reps);
    const max_iter = Number(raw.max_iter);

    const D_pi = asNumberArray(raw.D_pi) ?? [];
    const D_M  = asNumberArray(raw.D_M)  ?? [];

    if (!job_id || !Number.isFinite(thr) || !Number.isInteger(reps) || !Number.isInteger(max_iter)) {
      return NextResponse.json({ ok: false, error: "Invalid payload" }, { status: 400 });
    }

    const conn = db();
    await conn.execute(
      `INSERT INTO emtr_jobs
        (job_id, thr, reps, max_iter, d_pi_1, d_pi_2, d_pi_3, d_pi_4, d_m_1, d_m_2, d_m_3, d_m_4)
       VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
       ON DUPLICATE KEY UPDATE
         thr=VALUES(thr), reps=VALUES(reps), max_iter=VALUES(max_iter),
         d_pi_1=VALUES(d_pi_1), d_pi_2=VALUES(d_pi_2), d_pi_3=VALUES(d_pi_3), d_pi_4=VALUES(d_pi_4),
         d_m_1=VALUES(d_m_1),   d_m_2=VALUES(d_m_2),  d_m_3=VALUES(d_m_3),  d_m_4=VALUES(d_m_4)`,
      [
        job_id, thr, reps, max_iter,
        D_pi[0] ?? null, D_pi[1] ?? null, D_pi[2] ?? null, D_pi[3] ?? null,
        D_M[0]  ?? null, D_M[1]  ?? null, D_M[2]  ?? null, D_M[3]  ?? null,
      ]
    );

    return NextResponse.json({ ok: true });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
