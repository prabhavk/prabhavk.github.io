// app/api/wmw_comp_method/roots/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type Row = { root: string | null };

export async function GET(req: NextRequest) {
  const { searchParams } = new URL(req.url);
  const job = (searchParams.get("job") || "").trim();
  if (!job) {
    return NextResponse.json({ ok: false, error: "Missing ?job=<job_id>" }, { status: 400 });
  }

  try {
    const rows = await query<Row>(
      `
      SELECT DISTINCT r.root
        FROM emtr_init_final r
        JOIN emtr_jobs j ON j.job_id = r.job_id
       WHERE r.job_id = ?
         AND r.root IS NOT NULL
         AND r.root <> ''
         AND r.ll_final IS NOT NULL
         AND j.status = 'completed'
       ORDER BY r.root
      `,
      [job]
    );

    const roots = rows.map(r => r.root).filter((v): v is string => !!v);
    return NextResponse.json({ ok: true, roots });
  } catch (err) {
    const message = err instanceof Error ? err.message : String(err);
    return NextResponse.json({ ok: false, error: message }, { status: 500 });
  }
}
