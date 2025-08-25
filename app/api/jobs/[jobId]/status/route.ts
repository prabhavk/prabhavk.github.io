// app/api/jobs/[jobId]/status/route.ts
import { NextResponse } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Status = "started" | "completed" | "failed";

// --- small, typed helpers (no 'any') ---
function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isPromise<T = unknown>(x: unknown): x is Promise<T> {
  return isObject(x) && "then" in x && typeof (x as { then?: unknown }).then === "function";
}
async function getJobIdFromContext(ctx: unknown): Promise<string | null> {
  if (!isObject(ctx) || !("params" in ctx)) return null;
  const p = (ctx as { params: unknown }).params;

  if (isPromise<Record<string, unknown>>(p)) {
    const resolved = await p;
    return isObject(resolved) && typeof resolved.jobId === "string" ? resolved.jobId : null;
  }
  return isObject(p) && typeof p.jobId === "string" ? p.jobId : null;
}

export async function POST(req: Request, ctx: unknown) {
  try {
    const ct = req.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      return NextResponse.json(
        { ok: false, error: "Content-Type must be application/json" },
        { status: 415 }
      );
    }

    const jobId = await getJobIdFromContext(ctx);
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }

    const { status } = (await req.json().catch(() => ({}))) as { status?: Status };

    if (!status || !["started", "completed", "failed"].includes(status)) {
      return NextResponse.json({ ok: false, error: "Invalid status" }, { status: 400 });
    }

    const conn = db();
    await conn.execute(`UPDATE emtr_jobs SET status=? WHERE job_id=?`, [status, jobId]);

    return NextResponse.json({ ok: true });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
