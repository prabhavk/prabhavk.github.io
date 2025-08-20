// app/api/jobs/[jobId]/route.ts
import { NextRequest, NextResponse } from "next/server";
import { db } from "@/lib/pscale";

export const runtime = "edge";

type Status = "started" | "completed" | "failed";

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

export async function PATCH(
  req: NextRequest,
  { params }: { params: Promise<{ jobId: string }> } // Next 15: Promise
) {
  try {
    const { jobId } = await params;
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

    // Validate status if present
    if (status && !["started", "completed", "failed"].includes(status)) {
      return NextResponse.json({ ok: false, error: "Invalid status" }, { status: 400 });
    }

    // Build dynamic UPDATE
    const sets: string[] = [];
    const args: (string | number | null)[] = [];

    if (status) {
      sets.push("status = ?");
      args.push(status);
    }

    if (finished_at !== undefined) {
      sets.push("finished_at = ?");
      args.push(toMysqlTimestamp(finished_at)); // can be null to clear
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

export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ jobId: string }> }
) {
  try {
    const { jobId } = await params;
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
