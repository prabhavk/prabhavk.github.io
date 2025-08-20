// app/api/jobs/[jobId]/route.ts
import { NextRequest, NextResponse } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Status = "started" | "completed" | "failed";

export async function PATCH(
  req: NextRequest,
  { params }: { params: Promise<{ jobId: string }> } // ← Next 15: Promise
) {
  try {
    const { jobId } = await params; // ← await it
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }

    const body = await req.json().catch(() => null as unknown);
    const status = (body as { status?: Status })?.status;

    if (!status || !["started", "completed", "failed"].includes(status)) {
      return NextResponse.json({ ok: false, error: "Invalid status" }, { status: 400 });
    }

    const conn = db();
    await conn.execute(
      "UPDATE emtr_jobs SET status = ? WHERE job_id = ?",
      [status, jobId]
    );

    return NextResponse.json({ ok: true, jobId, status });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}

export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ jobId: string }> } // ← Next 15: Promise
) {
  try {
    const { jobId } = await params; // ← await it
    if (!jobId) {
      return NextResponse.json({ ok: false, error: "Missing jobId" }, { status: 400 });
    }
    const conn = db();
    const res = await conn.execute(
      "SELECT * FROM emtr_jobs WHERE job_id = ? LIMIT 1",
      [jobId]
    );
    return NextResponse.json({ ok: true, job: (res.rows as Record<string, unknown>[])[0] ?? null });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
