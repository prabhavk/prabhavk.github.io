// app/api/jobs/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Status = "started" | "completed" | "failed";

interface JobPayload {
  job_id?: string;
  thr?: number | string;
  reps?: number | string;
  max_iter?: number | string;
  D_pi?: Array<number | string>; // length 4
  D_M?: Array<number | string>;  // length 4
  status?: Status;
  secret?: string;               // optional body fallback
}

const DEFAULT_D_PI = [100, 100, 100, 100] as const;
const DEFAULT_D_M  = [100, 2, 2, 2] as const;

function toNumber(x: unknown): number {
  const n = typeof x === "string" ? Number(x.trim()) : Number(x);
  return Number.isFinite(n) ? n : NaN;
}

function asNumberArray(x: unknown): number[] | null {
  if (!Array.isArray(x)) return null;
  const arr = x.map((v) => toNumber(v));
  return arr.every((v) => Number.isFinite(v)) ? arr : null;
}

function isValidDirichlet4(a: number[] | null): a is number[] {
  return !!a && a.length === 4 && a.every((v) => Number.isFinite(v) && v > 0);
}

function isAllowedBySecret(secret?: string | null): boolean {
  const raw = process.env.EMTR_ALLOWED_SECRETS;
  if (!raw) return true;
  const allowed = raw.split(/[,\s]+/g).map(s => s.trim()).filter(Boolean);
  if (allowed.length === 0) return true;
  return !!secret && allowed.includes(secret);
}

export async function POST(req: NextRequest) {
  try {
    const ct = (req.headers.get("content-type") || "").toLowerCase();
    if (!ct.includes("application/json")) {
      return NextResponse.json({ ok: false, error: "Content-Type must be application/json" }, { status: 415 });
    }

    const raw = (await req.json()) as Partial<JobPayload>;

    // single source of truth for job_id
    let job_id = String(raw.job_id ?? "").trim();
    if (!job_id) job_id = `job-${Date.now()}`;

    const thr = toNumber(raw.thr);
    const reps = toNumber(raw.reps);
    const max_iter = toNumber(raw.max_iter);

    if (!Number.isFinite(thr) || thr <= 0) {
      return NextResponse.json({ ok: false, error: "thr must be a positive number" }, { status: 400 });
    }
    if (!Number.isInteger(reps) || reps < 1) {
      return NextResponse.json({ ok: false, error: "reps must be a positive integer" }, { status: 400 });
    }
    if (!Number.isInteger(max_iter) || max_iter < 1) {
      return NextResponse.json({ ok: false, error: "max_iter must be a positive integer" }, { status: 400 });
    }

    const D_pi_in = asNumberArray(raw.D_pi);
    const D_M_in  = asNumberArray(raw.D_M);
    const D_pi = isValidDirichlet4(D_pi_in) ? D_pi_in : [...DEFAULT_D_PI];
    const D_M  = isValidDirichlet4(D_M_in)  ? D_M_in  : [...DEFAULT_D_M];

    const status: Status =
      raw.status === "completed" || raw.status === "failed" || raw.status === "started"
        ? raw.status
        : "started";

    // accept secret via header OR body
    const suppliedSecret = req.headers.get("x-emtr-secret") ?? raw.secret ?? null;
    if (!isAllowedBySecret(suppliedSecret)) {
      return NextResponse.json({ ok: false, error: "Forbidden: invalid or missing secret" }, { status: 403 });
    }

    const conn = db();
    await conn.execute(
      `INSERT INTO emtr_jobs
         (job_id, status, thr, reps, max_iter,
          d_pi_1, d_pi_2, d_pi_3, d_pi_4,
          d_m_1,  d_m_2,  d_m_3,  d_m_4)
       VALUES (?,?,?,?,?,
               ?,?,?,?,
               ?,?,?,?)
       ON DUPLICATE KEY UPDATE
         status=VALUES(status),
         thr=VALUES(thr),
         reps=VALUES(reps),
         max_iter=VALUES(max_iter),
         d_pi_1=VALUES(d_pi_1), d_pi_2=VALUES(d_pi_2), d_pi_3=VALUES(d_pi_3), d_pi_4=VALUES(d_pi_4),
         d_m_1=VALUES(d_m_1),   d_m_2=VALUES(d_m_2),  d_m_3=VALUES(d_m_3),  d_m_4=VALUES(d_m_4)`,
      [
        job_id, status, thr, reps, max_iter,
        D_pi[0], D_pi[1], D_pi[2], D_pi[3],
        D_M[0],  D_M[1],  D_M[2],  D_M[3],
      ]
    );

    // return shape that matches your client: data.job_id
    return NextResponse.json({
      ok: true,
      job_id,
      status,
      thr,
      reps,
      max_iter,
      D_pi,
      D_M,
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
