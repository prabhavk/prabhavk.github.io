// app/api/familyjoining/route.ts
import { NextRequest, NextResponse } from "next/server";
import { query } from "@/lib/db";

export const runtime = "nodejs";

type EdgeTuple = [string, string, number];

function toNumber(x: unknown): number | null {
  const n = Number(x);
  return Number.isFinite(n) ? n : null;
}
function isRecord(v: unknown): v is Record<string, unknown> {
  return typeof v === "object" && v !== null;
}
function safeParseJSON(s: string): unknown | null {
  try { return JSON.parse(s); } catch { return null; }
}

function coerceEdges(body: Record<string, unknown>): EdgeTuple[] | null {
  const candidate =
    body["edge_list"] ??
    body["edges"] ??
    body["edgeObjects"] ??
    body["edge_objects"] ??
    null;

  let raw: unknown = candidate;

  if (typeof raw === "string") {
    const parsed = safeParseJSON(raw);
    if (parsed === null) return null;
    raw = parsed;
  }

  if (!Array.isArray(raw)) return null;

  const out: EdgeTuple[] = [];
  for (const e of raw as unknown[]) {
    if (Array.isArray(e) && e.length >= 3) {
      const u = typeof e[0] === "string" ? e[0] : String(e[0] ?? "");
      const v = typeof e[1] === "string" ? e[1] : String(e[1] ?? "");
      const w = Number(e[2]);
      if (!u || !v || !Number.isFinite(w)) return null;
      out.push([u, v, w]);
      continue;
    }
    if (isRecord(e)) {
      const u = typeof e["u"] === "string" ? e["u"] as string : String(e["u"] ?? "");
      const v = typeof e["v"] === "string" ? e["v"] as string : String(e["v"] ?? "");
      const w = Number(e["w"]);
      if (!u || !v || !Number.isFinite(w)) return null;
      out.push([u, v, w]);
      continue;
    }
    return null;
  }
  return out.length ? out : null;
}

function countLeaves(edges: EdgeTuple[]): number {
  const names = new Set<string>();
  for (const [u, v] of edges) { names.add(u); names.add(v); }
  return names.size;
}

export async function POST(req: NextRequest) {
  try {
    const headers = req.headers;
    const bodyText = await req.text().catch(() => "");
    const parsed: unknown = bodyText ? safeParseJSON(bodyText) : {};
    const b: Record<string, unknown> = isRecord(parsed) ? parsed : {};

    // job_id
    const jobId =
      headers.get("x-job-id") ??
      (typeof b["job_id"] === "string" ? (b["job_id"] as string) : null);

    if (!jobId) {
      return NextResponse.json(
        { ok: false, error: "Missing job_id" },
        { status: 400 }
      );
    }

    // edges
    const edges = coerceEdges(b);
    if (!edges) {
      return NextResponse.json(
        { ok: false, error: "Missing or invalid edge list" },
        { status: 400 }
      );
    }

    // fields (rep is parsed but not stored)
    const method  = typeof b["method"] === "string" ? (b["method"] as string) : "na";
    const label   = typeof b["label"]  === "string" ? (b["label"]  as string) : "final";
    const source  = typeof b["source"] === "string" ? (b["source"] as string) : "familyjoining";
    const rep     = toNumber(b["rep"]);      // informational only
    const iter    = toNumber(b["iter"]);
    const epsilon = toNumber(b["epsilon"]);
    const root    = typeof b["root"] === "string" ? (b["root"] as string) : "";

    const n_edges  = Number.isFinite(Number(b["n_edges"]))  ? Number(b["n_edges"])  : edges.length;
    const n_leaves = Number.isFinite(Number(b["n_leaves"])) ? Number(b["n_leaves"]) : countLeaves(edges);

    // Upsert header (keyed by job_id, method, label)
    await query(
      `INSERT INTO emtr_fj_runs
         (job_id, method, label, source, iter, epsilon, root, n_edges, n_leaves, raw_json)
       VALUES
         (?, ?, ?, ?, ?, ?, ?, ?, ?, CAST(? AS JSON))
       ON DUPLICATE KEY UPDATE
         source   = VALUES(source),
         iter     = VALUES(iter),
         epsilon  = VALUES(epsilon),
         root     = VALUES(root),
         n_edges  = VALUES(n_edges),
         n_leaves = VALUES(n_leaves),
         raw_json = VALUES(raw_json)`,
      [
        jobId,
        method,
        label,
        source,
        iter,
        epsilon,
        root,
        n_edges,
        n_leaves,
        JSON.stringify(b),
      ]
    );

    // Replace edges for (job_id, method, label)
    await query(
      `DELETE FROM emtr_fj_edges
       WHERE job_id = ? AND method = ? AND label = ?`,
      [jobId, method, label]
    );

    if (edges.length) {
      const chunkSize = 1000;
      for (let i = 0; i < edges.length; i += chunkSize) {
        const chunk = edges.slice(i, i + chunkSize);
        const placeholders = chunk.map(() => "(?, ?, ?, ?, ?, ?)").join(", ");
        const params: (string | number)[] = [];
        for (const [u, v, w] of chunk) {
          params.push(jobId, method, label, u, v, w);
        }
        await query(
          `INSERT INTO emtr_fj_edges
             (job_id, method, label, u, v, w)
           VALUES ${placeholders}
           ON DUPLICATE KEY UPDATE
             w = VALUES(w)`,
          params
        );
      }
    }

    return NextResponse.json({
      ok: true,
      data: {
        job_id: jobId,
        source, method, label, rep, iter, epsilon, root,
        n_edges, n_leaves,
        edges,
      },
    });
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : String(e);
    return NextResponse.json(
      { ok: false, error: `DB error: ${msg}` },
      { status: 500 }
    );
  }
}
