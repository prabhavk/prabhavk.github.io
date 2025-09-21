// app/api/ui-inputs/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db"; // or use your { db } helper if you prefer Edge

export const runtime = "nodejs"; // or "edge" if your db() helper supports it

type Payload = {
  job_id: string;
  settings_version: number;
  use_default_seqs: boolean;
  dna_filename: string | null;
  thr: number;
  reps: number;
  max_iter: number;
  require_alpha_ge1: boolean;
  include_parsimony: boolean;
  include_hss: boolean;
  D_pi: [number, number, number, number];
  D_M: [number, number, number, number];
  pattern_counts: { coarse: number; medium: number; fine: number };
  raw_json: unknown; // full payload you send to the worker
};

function isNum(x: unknown): x is number {
  return typeof x === "number" && Number.isFinite(x);
}
function isBool(x: unknown): x is boolean {
  return typeof x === "boolean";
}
function isVec4(a: unknown): a is [number, number, number, number] {
  return Array.isArray(a) && a.length === 4 && a.every(isNum);
}

export async function POST(req: NextRequest) {
  let b: Payload;
  try {
    b = await req.json();
  } catch {
    return NextResponse.json({ ok: false, error: "Invalid JSON" }, { status: 400 });
  }
  if (!b?.job_id) {
    return NextResponse.json({ ok: false, error: "job_id required" }, { status: 400 });
  }
  if (
    !isNum(b.settings_version) ||
    !isNum(b.thr) ||
    !isNum(b.reps) ||
    !isNum(b.max_iter) ||
    !isBool(b.require_alpha_ge1) ||
    !isBool(b.include_parsimony) ||
    !isBool(b.include_hss) ||
    !isVec4(b.D_pi) ||
    !isVec4(b.D_M) ||
    !isNum(b.pattern_counts?.coarse) ||
    !isNum(b.pattern_counts?.medium) ||
    !isNum(b.pattern_counts?.fine)
  ) {
    return NextResponse.json({ ok: false, error: "Missing/invalid fields" }, { status: 400 });
  }

  // Upsert (PlanetScale supports INSERT ... ON DUPLICATE KEY UPDATE)
  const sql = `
    INSERT INTO emtr_ui_inputs
      (job_id, settings_version, use_default_seqs, dna_filename,
       thr, reps, max_iter, require_alpha_ge1, include_parsimony, include_hss,
       D_pi, D_M, pattern_counts, raw_json)
    VALUES
      (?, ?, ?, ?,
       ?, ?, ?, ?, ?, ?,
       CAST(? AS JSON), CAST(? AS JSON), CAST(? AS JSON), CAST(? AS JSON))
    ON DUPLICATE KEY UPDATE
       settings_version=VALUES(settings_version),
       use_default_seqs=VALUES(use_default_seqs),
       dna_filename=VALUES(dna_filename),
       thr=VALUES(thr),
       reps=VALUES(reps),
       max_iter=VALUES(max_iter),
       require_alpha_ge1=VALUES(require_alpha_ge1),
       include_parsimony=VALUES(include_parsimony),
       include_hss=VALUES(include_hss),
       D_pi=VALUES(D_pi),
       D_M=VALUES(D_M),
       pattern_counts=VALUES(pattern_counts),
       raw_json=VALUES(raw_json)
  `;
  await query(sql, [
    b.job_id,
    b.settings_version,
    b.use_default_seqs ? 1 : 0,
    b.dna_filename,
    b.thr,
    b.reps,
    b.max_iter,
    b.require_alpha_ge1 ? 1 : 0,
    b.include_parsimony ? 1 : 0,
    b.include_hss ? 1 : 0,
    JSON.stringify(b.D_pi),
    JSON.stringify(b.D_M),
    JSON.stringify(b.pattern_counts),
    JSON.stringify(b.raw_json),
  ]);

  return NextResponse.json({ ok: true });
}
