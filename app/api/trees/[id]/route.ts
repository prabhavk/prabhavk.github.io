// app/api/trees/[id]/route.ts
import { NextResponse } from "next/server";
import { q } from "@/lib/db";

type YesNo = 0 | 1;

// Satisfies q<T extends Record<string, unknown>>
interface EmtrTreeRow {
  [key: string]: unknown;
  id: number;
  job_id: string;
  label: string | null;
  method: string;
  source: string;
  epsilon: number | null;
  tree: string;       // newick
  created_at: string; // TIMESTAMP usually as string
  is_current: YesNo;
}

type JobIdRow = { [k: string]: unknown; job_id: string };

function parseIdFromUrl(req: Request): number | null {
  try {
    const { pathname } = new URL(req.url);
    // pathname like: /api/trees/123
    const segs = pathname.split("/").filter(Boolean);
    const idStr = segs[segs.length - 1] ?? "";
    const n = Number.parseInt(idStr, 10);
    return Number.isFinite(n) ? n : null;
  } catch {
    return null;
  }
}

export async function GET(req: Request) {
  const id = parseIdFromUrl(req);
  if (id === null) return NextResponse.json({ error: "invalid id" }, { status: 400 });

  const rows = await q<EmtrTreeRow>(`SELECT * FROM emtr_trees WHERE id=?`, [id]);
  const row = rows[0];
  if (!row) return NextResponse.json({ error: "not found" }, { status: 404 });
  return NextResponse.json(row);
}

export async function PATCH(req: Request) {
  const id = parseIdFromUrl(req);
  if (id === null) return NextResponse.json({ error: "invalid id" }, { status: 400 });

  const body = (await req.json().catch(() => ({}))) as Record<string, unknown>;
  const is_current = Number(body.is_current ?? 0) | 0;

  if (is_current === 1) {
    const rows = await q<JobIdRow>(`SELECT job_id FROM emtr_trees WHERE id=?`, [id]);
    const found = rows[0];
    if (!found) return NextResponse.json({ error: "not found" }, { status: 404 });

    await q(`UPDATE emtr_trees SET is_current=0 WHERE job_id=?`, [found.job_id as string]);
    await q(`UPDATE emtr_trees SET is_current=1 WHERE id=?`, [id]);
    return NextResponse.json({ ok: true });
  }

  // Extend with other editable fields if needed
  return NextResponse.json({ ok: true });
}

export async function DELETE(req: Request) {
  const id = parseIdFromUrl(req);
  if (id === null) return NextResponse.json({ error: "invalid id" }, { status: 400 });

  await q(`DELETE FROM emtr_trees WHERE id=?`, [id]);
  return NextResponse.json({ ok: true });
}
