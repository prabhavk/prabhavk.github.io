// app/api/jobs/[jobId]/status/route.ts
import { NextResponse, type NextRequest } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Status = "started" | "completed" | "failed";

export async function POST(
  req: NextRequest,
  { params }: { params: Promise<{ jobId: string }> }  // <- Promise
) {
  try {
    const ct = req.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      return NextResponse.json({ ok: false, error: "Content-Type must be application/json" }, { status: 415 });
    }

    const { jobId } = await params;  // <- await it
    const { status } = (await req.json()) as { status?: Status };

    if (!status || !["started", "completed", "failed"].includes(status)) {
      return NextResponse.json({ ok: false, error: "Invalid status" }, { status: 400 });
    }

    const conn = db();
    await conn.execute(`UPDATE emtr_jobs SET status=? WHERE job_id=?`, [status, jobId]);

    return NextResponse.json({ ok: true });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ ok: false, error: msg }, { status: 500 });
  }
}
