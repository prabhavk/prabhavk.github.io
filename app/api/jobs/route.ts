// app/api/jobs/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Status = "started" | "completed" | "failed";

interface JobPayload {
  job_id: string;
  thr: number;
  reps: number;
  max_iter: number;
  D_pi?: number[]; // length 4
  D_M?: number[];  // length 4
  status?: Status; // defaults to "started"
}

function asNumberArray(x: unknown): number[] | null {
  if (!Array.isArray(x)) return null;
  const arr = x.map((v) => Number(v));
  return arr.every((v) => Number.isFinite(v)) ? arr : null;
}

function isValidDirichlet4(a: number[] | null): a is number[] {
  return !!a && a.length === 4 && a.every((v) => Number.isFinite(v) && v > 0);
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

    const D_pi = asNumberArray(raw.D_pi);
    const D_M  = asNumberArray(raw.D_M);

    const status: Status =
      raw.status === "completed" || raw.status === "failed" || raw.status === "started"
        ? raw.status
        : "started";

    // Basic numeric checks
    if (!job_id || !Number.isFinite(thr) || thr <= 0) {
      return NextResponse.json({ ok: false, error: "Invalid job_id or thr" }, { status: 400 });
    }
    if (!Number.isInteger(reps) || reps < 1) {
      return NextResponse.json({ ok: false, error: "reps must be a positive integer" }, { status: 400 });
    }
    if (!Number.isInteger(max_iter) || max_iter < 1) {
      return NextResponse.json({ ok: false, error: "max_iter must be a positive integer" }, { status: 400 });
    }

    // Dirichlet priors (optional, but if provided must be length-4 and > 0)
    if (raw.D_pi && !isValidDirichlet4(D_pi)) {
      return NextResponse.json({ ok: false, error: "D_pi must be an array of 4 positive numbers" }, { status: 400 });
    }
    if (raw.D_M && !isValidDirichlet4(D_M)) {
      return NextResponse.json({ ok: false, error: "D_M must be an array of 4 positive numbers" }, { status: 400 });
    }

    const conn = db();
    await conn.execute(
      `INSERT INTO emtr_jobs
        (job_id, status, thr, reps, max_iter,
         d_pi_1, d_pi_2, d_pi_3, d_pi_4,
         d_m_1,  d_m_2,  d_m_3,  d_m_4)
       VALUES (?,?,?,?, ?,?,?,?, ?,?,?,?, ?)
       ON DUPLICATE KEY UPDATE
         status=VALUES(status),
         thr=VALUES(thr),
         reps=VALUES(reps),
         max_iter=VALUES(max_iter),
         d_pi_1=VALUES(d_pi_1), d_pi_2=VALUES(d_pi_2), d_pi_3=VALUES(d_pi_3), d_pi_4=VALUES(d_pi_4),
         d_m_1=VALUES(d_m_1),   d_m_2=VALUES(d_m_2),  d_m_3=VALUES(d_m_3),  d_m_4=VALUES(d_m_4)`,
      [
        job_id, status, thr, reps, max_iter,
        D_pi?.[0] ?? null, D_pi?.[1] ?? null, D_pi?.[2] ?? null, D_pi?.[3] ?? null,
        D_M?.[0]  ?? null, D_M?.[1]  ?? null, D_M?.[2]  ?? null, D_M?.[3]  ?? null,
      ]
    );

    return NextResponse.json({
      ok: true,
      job: {
        job_id, status, thr, reps, max_iter,
        D_pi: D_pi ?? null,
        D_M:  D_M  ?? null,
      },
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
