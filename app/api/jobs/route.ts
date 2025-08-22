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

const DEFAULT_D_PI = [100, 100, 100, 100] as const;
const DEFAULT_D_M  = [100, 2, 2, 2] as const;

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

    const D_pi_in = asNumberArray(raw.D_pi);
    const D_M_in  = asNumberArray(raw.D_M);

    const status: Status =
      raw.status === "completed" || raw.status === "failed" || raw.status === "started"
        ? raw.status
        : "started";

    if (!job_id || job_id.length > 128) {
      return NextResponse.json({ ok: false, error: "Invalid job_id" }, { status: 400 });
    }
    if (!Number.isFinite(thr) || thr <= 0) {
      return NextResponse.json({ ok: false, error: "thr must be a positive number" }, { status: 400 });
    }
    if (!Number.isInteger(reps) || reps < 1) {
      return NextResponse.json({ ok: false, error: "reps must be a positive integer" }, { status: 400 });
    }
    if (!Number.isInteger(max_iter) || max_iter < 1) {
      return NextResponse.json({ ok: false, error: "max_iter must be a positive integer" }, { status: 400 });
    }

    if (raw.D_pi && !isValidDirichlet4(D_pi_in)) {
      return NextResponse.json(
        { ok: false, error: "D_pi must be an array of 4 positive numbers" },
        { status: 400 }
      );
    }
    if (raw.D_M && !isValidDirichlet4(D_M_in)) {
      return NextResponse.json(
        { ok: false, error: "D_M must be an array of 4 positive numbers" },
        { status: 400 }
      );
    }
    const D_pi = isValidDirichlet4(D_pi_in) ? D_pi_in : [...DEFAULT_D_PI];
    const D_M  = isValidDirichlet4(D_M_in)  ? D_M_in  : [...DEFAULT_D_M];

    const conn = db();
    await conn.execute(
      `INSERT INTO emtr_jobs
         (job_id, status, thr, reps, max_iter, d_pi, d_m)
       VALUES (?,?,?,?,?, CAST(? AS JSON), CAST(? AS JSON))
       ON DUPLICATE KEY UPDATE
         status=VALUES(status),
         thr=VALUES(thr),
         reps=VALUES(reps),
         max_iter=VALUES(max_iter),
         d_pi=VALUES(d_pi),
         d_m=VALUES(d_m)`,
      [
        job_id, status, thr, reps, max_iter,
        JSON.stringify(D_pi),
        JSON.stringify(D_M),
      ]
    );

    return NextResponse.json({
      ok: true,
      job: { job_id, status, thr, reps, max_iter, D_pi, D_M },
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
