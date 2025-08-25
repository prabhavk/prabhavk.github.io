// app/api/jobs/[jobId]/route.ts
import { NextResponse } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Status = "started" | "completed" | "failed";

// ---------- small helpers (no 'any' annotations) ----------
function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isPromise<T = unknown>(x: unknown): x is Promise<T> {
  // Safe check without annotating 'any' at the signature level
  return isObject(x) && "then" in x && typeof (x as { then?: unknown }).then === "function";
}

async function getJobIdFromContext(ctx: unknown): Promise<string | null> {
  if (!isObject(ctx) || !("params" in ctx)) return null;
  const p = (ctx as { params: unknown }).params;

  if (isPromise<Record<string, unknown>>(p)) {
    const resolved = await p;
    const jid = isObject(resolved) && typeof resolved.jobId === "string" ? resolved.jobId : null;
    return jid;
  }

  const jid = isObject(p) && typeof p.jobId === "string" ? p.jobId : null;
  return jid;
}

function toMysqlTimestamp(x: unknown): string | null {
  if (x == null) return null;
  const d =
    typeof x === "number"
      ? new Date(x)
      : typeof x === "string"
      ? new Date(x)
      : null;
  if (!d || Number.isNaN(d.getTime())) return null;
  // "YYYY-MM-DD HH:MM:SS"
  return new Date(d.getTime() - d.getTimezoneOffset() * 60000)
    .toISOString()
    .slice(0, 19)
    .replace("T", " ");
}

// ---------- PATCH ----------
export async function PATCH(req: Request, ctx: unknown) {
  try {
    const jobId = await getJobIdFromContext(ctx);
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }

    const body = (await req.json().catch(() => null)) as
      | {
          status?: Status;
          finished_at?: string | number | null;
          dur_minutes?: number | null;
        }
      | null;

    if (!body) {
      return NextResponse.json({ ok: false, error: "Invalid JSON body" }, { status: 400 });
    }

    const { status, finished_at, dur_minutes } = body;

    if (status && !["started", "completed", "failed"].includes(status)) {
      return NextResponse.json({ ok: false, error: "Invalid status" }, { status: 400 });
    }

    const sets: string[] = [];
    const args: (string | number | null)[] = [];

    if (status) {
      sets.push("status = ?");
      args.push(status);
    }
    if (finished_at !== undefined) {
      sets.push("finished_at = ?");
      args.push(toMysqlTimestamp(finished_at)); // null clears column
    }
    if (dur_minutes !== undefined) {
      const dur =
        dur_minutes === null
          ? null
          : Number.isFinite(Number(dur_minutes)) && Number(dur_minutes) >= 0
          ? Math.floor(Number(dur_minutes))
          : null;
      sets.push("dur_minutes = ?");
      args.push(dur);
    }

    if (sets.length === 0) {
      return NextResponse.json(
        { ok: false, error: "Nothing to update; provide at least one of status, finished_at, dur_minutes" },
        { status: 400 }
      );
    }

    const sql = `UPDATE emtr_jobs SET ${sets.join(", ")} WHERE job_id = ?`;
    args.push(jobId);

    const conn = db();
    await conn.execute(sql, args);

    return NextResponse.json({
      ok: true,
      jobId,
      updated: {
        ...(status ? { status } : {}),
        ...(finished_at !== undefined ? { finished_at } : {}),
        ...(dur_minutes !== undefined ? { dur_minutes } : {}),
      },
    });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}

// ---------- GET ----------
export async function GET(_req: Request, ctx: unknown) {
  try {
    const jobId = await getJobIdFromContext(ctx);
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }

    const conn = db();
    const res = await conn.execute("SELECT * FROM emtr_jobs WHERE job_id = ? LIMIT 1", [jobId]);
    return NextResponse.json({ ok: true, job: (res.rows as Record<string, unknown>[])[0] ?? null });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
