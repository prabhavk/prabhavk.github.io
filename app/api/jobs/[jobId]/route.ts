// app/api/jobs/[jobId]/route.ts
import { NextResponse } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Status = "started" | "completed" | "failed" | "blocked";

/* ------------------------ small helpers ------------------------ */
function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}

function isPromise<T = unknown>(x: unknown): x is Promise<T> {
  return isObject(x) && "then" in x && typeof (x as { then?: unknown }).then === "function";
}

/** Next.js can pass params or a promise of params on edge routes */
async function getJobIdFromContext(ctx: unknown): Promise<string | null> {
  if (!isObject(ctx) || !("params" in ctx)) return null;
  const p = (ctx as { params: unknown }).params;

  if (isPromise<Record<string, unknown>>(p)) {
    const resolved = await p;
    const jid = isObject(resolved) && typeof resolved.jobId === "string" ? resolved.jobId : null;
    return jid;
  }
  const jid = isObject(p) && typeof (p as Record<string, unknown>).jobId === "string"
    ? ((p as Record<string, string>).jobId)
    : null;
  return jid;
}

/** Format as naive UTC DATETIME 'YYYY-MM-DD HH:MM:SS' for MySQL/PlanetScale */
function toMysqlTimestampUTC(x: unknown): string | null {
  if (x == null) return null;
  const d = typeof x === "number" ? new Date(x)
        : typeof x === "string" ? new Date(x)
        : null;
  if (!d || Number.isNaN(d.getTime())) return null;
  return d.toISOString().slice(0, 19).replace("T", " ");
}

/** Optional gating via EMTR_ALLOWED_SECRETS (comma/space separated) */
function isAllowedBySecret(secret: string | null): boolean {
  const raw = process.env.EMTR_ALLOWED_SECRETS;
  if (!raw) return true;
  const allow = raw.split(/[,\s]+/g).map(s => s.trim()).filter(Boolean);
  if (allow.length === 0) return true;
  return !!secret && allow.includes(secret);
}

function getSuppliedSecretFromHeaders(req: Request): string | null {
  const hdr = req.headers.get("x-emtr-secret");
  const auth = req.headers.get("authorization") || "";
  const fromAuth = auth.toLowerCase().startsWith("bearer ") ? auth.slice(7).trim() : "";
  return (hdr && hdr.trim()) || (fromAuth && fromAuth.trim()) || null;
}

/* ------------------------------- PATCH ------------------------------- */
export async function PATCH(req: Request, ctx: unknown) {
  try {
    const job_id = await getJobIdFromContext(ctx);
    if (!job_id) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }

    // Optional secret gating (match /api/jobs behavior)
    if (!isAllowedBySecret(getSuppliedSecretFromHeaders(req))) {
      return NextResponse.json({ ok: false, error: "Forbidden: invalid or missing secret" }, { status: 403 });
    }

    const body = (await req.json().catch(() => null)) as
      | {
          status?: Status;
          started_at?: string | number | null;
          finished_at?: string | number | null;
          dur_minutes?: number | null;
          blocked_reason?: string | null;
        }
      | null;

    if (!body) {
      return NextResponse.json({ ok: false, error: "Invalid JSON body" }, { status: 400 });
    }

    const { status, started_at, finished_at, dur_minutes, blocked_reason } = body;

    if (status && !["started", "completed", "failed", "blocked"].includes(status)) {
      return NextResponse.json({ ok: false, error: "Invalid status" }, { status: 400 });
    }

    const sets: string[] = [];
    const args: (string | number | null)[] = [];

    if (status) {
      sets.push("status = ?");
      args.push(status);
    }

    if (started_at !== undefined) {
      const startedAtSql = toMysqlTimestampUTC(started_at); // null clears
      sets.push("started_at = ?");
      args.push(startedAtSql);
    }

    let finishedAtSql: string | null | undefined = undefined;
    if (finished_at !== undefined) {
      finishedAtSql = toMysqlTimestampUTC(finished_at); // null clears column when null provided
      sets.push("finished_at = ?");
      args.push(finishedAtSql);
    }

    let userSetDuration = false;
    if (dur_minutes !== undefined) {
      userSetDuration = true;
      const dur =
        dur_minutes === null
          ? null
          : Number.isFinite(Number(dur_minutes)) && Number(dur_minutes) >= 0
          ? Number(dur_minutes)
          : null;
      sets.push("dur_minutes = ?");
      args.push(dur);
    }

    if (blocked_reason !== undefined) {
      // Store as-is or trim if you prefer to enforce a max length
      const reason = blocked_reason == null ? null : String(blocked_reason).slice(0, 2048);
      sets.push("blocked_reason = ?");
      args.push(reason);
    }

    // If finished_at provided and dur_minutes not provided, compute automatically in SQL
    if (finished_at !== undefined && !userSetDuration) {
      sets.push("dur_minutes = ROUND(TIMESTAMPDIFF(SECOND, created_at, finished_at) / 60.0, 2)");
      // (no bound arg)
    }

    if (sets.length === 0) {
      return NextResponse.json(
        { ok: false, error: "Nothing to update; provide one of status, started_at, finished_at, dur_minutes, blocked_reason" },
        { status: 400 }
      );
    }

    const sql = `UPDATE emtr_jobs SET ${sets.join(", ")} WHERE job_id = ?`;
    args.push(job_id);

    const conn = db();
    await conn.execute(sql, args);

    return NextResponse.json({
      ok: true,
      job_id,
      updated: {
        ...(status !== undefined ? { status } : {}),
        ...(started_at !== undefined ? { started_at } : {}),
        ...(finished_at !== undefined ? { finished_at } : {}),
        ...(dur_minutes !== undefined ? { dur_minutes } : {}),
        ...(blocked_reason !== undefined ? { blocked_reason } : {}),
      },
    });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}

/* -------------------------------- GET -------------------------------- */
export async function GET(_req: Request, ctx: unknown) {
  try {
    const job_id = await getJobIdFromContext(ctx);
    if (!job_id) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }

    const conn = db();
    const res = await conn.execute("SELECT * FROM emtr_jobs WHERE job_id = ? LIMIT 1", [job_id]);
    const row = (res.rows as Record<string, unknown>[])[0] ?? null;

    return NextResponse.json({ ok: true, job: row });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
