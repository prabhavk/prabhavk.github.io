// app/api/jobs/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

export async function POST(req: NextRequest) {
  try {
    const ct = req.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      return NextResponse.json({ ok: false, error: "Content-Type must be application/json" }, { status: 415 });
    }
    const body = await req.json();
    const job_id = String(body.job_id ?? "").trim();
    const thr = Number(body.thr);
    const reps = Number(body.reps);
    const max_iter = Number(body.max_iter);
    const D_pi = Array.isArray(body.D_pi) ? body.D_pi.map(Number) : [];
    const D_M  = Array.isArray(body.D_M)  ? body.D_M.map(Number)  : [];

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
  } catch (e) {
    const msg = e && typeof e === "object" && "message" in e ? String((e as any).message) : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
