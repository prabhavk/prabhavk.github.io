// app/api/trees/route.ts
import { NextResponse } from "next/server";
import { query as q } from "@/lib/db";
import crypto from "node:crypto";

export const runtime = "nodejs"; // needed for node:crypto

type YesNo = 0 | 1;

type Edge = [string, string, number];

type TreeRow = {
  [key: string]: unknown;
  id: number;
  job_id: string;
  source: string;
  method: string;
  label: string | null;
  format: "newick" | "nexus" | "json" | string;
  tree?: string | null;           // present when full=1
  tree_labeled?: string | null;   // present when full=1
  edge_list?: Edge[] | null;      // present when full=1 (driver may return as JSON string)
  tree_preview: string;           // always present
  tree_len: number | null;
  rep: number | null;
  iter: number | null;
  epsilon: number | null;
  n_leaves: number | null;
  n_edges: number | null;
  root: string | null;
  is_current: YesNo;
  created_at: string;
  updated_at: string;
};

function truthy(v: string | null): boolean {
  return v !== null && /^(1|true|yes|on)$/i.test(v);
}
function clampInt(n: unknown, def: number, min: number, max: number): number {
  const x = Number.parseInt(String(n ?? ""), 10);
  if (!Number.isFinite(x)) return def;
  return Math.max(min, Math.min(x, max));
}
function coerceYesNo(v: unknown): YesNo {
  if (typeof v === "number") return v ? 1 : 0;
  if (typeof v === "boolean") return v ? 1 : 0;
  if (typeof v === "string") return truthy(v) ? 1 : 0;
  return 0;
}
function toNullableNumber(v: unknown): number | null {
  if (v === null || v === undefined || v === "") return null;
  const n = Number(v);
  return Number.isFinite(n) ? n : null;
}

function normalizeEdgeList(raw: unknown): Edge[] | null {
  let x = raw as unknown;
  if (typeof x === "string") {
    try { x = JSON.parse(x); } catch { return null; }
  }
  if (!Array.isArray(x)) return null;
  const out: Edge[] = [];
  for (const e of x) {
    if (!Array.isArray(e) || e.length !== 3) return null;
    const a = e[0], b = e[1], w = Number(e[2]);
    if (typeof a !== "string" || typeof b !== "string" || !Number.isFinite(w)) return null;
    out.push([a, b, w]);
  }
  return out;
}

function inferLeafCountFromEdges(edges: Edge[]): number {
  // Per your rule: leaf nodes do not have h_ in their name; only internal nodes do.
  const names = new Set<string>();
  for (const [u, v] of edges) { names.add(u); names.add(v); }
  let cnt = 0;
  for (const name of names) {
    if (!/^h_/.test(name)) cnt++;
  }
  return cnt;
}

/* ---------------- GET ---------------- */

export async function GET(req: Request) {
  try {
    const url = new URL(req.url);

    const job_id = (url.searchParams.get("job_id") || url.searchParams.get("job") || "").trim();
    const currentOnly = truthy(url.searchParams.get("current"));
    const full = truthy(url.searchParams.get("full"));
    const fullInt = full ? 1 : 0;

    const limit = clampInt(url.searchParams.get("limit"), 50, 1, 200);
    const offset = clampInt(url.searchParams.get("offset"), 0, 0, 1_000_000);

    if (!job_id) {
      return NextResponse.json({ error: "job_id required" }, { status: 400 });
    }

    const where = `WHERE job_id = ? ${currentOnly ? "AND is_current=1" : ""}`;
    const rows = await q<TreeRow>(
      `
      SELECT
        id, job_id, source, method, label, format,
        IF(?=1, emtr_trees.tree, NULL)          AS tree,
        IF(?=1, emtr_trees.tree_labeled, NULL)  AS tree_labeled,
        IF(?=1, emtr_trees.edge_list, NULL)     AS edge_list,
        LEFT(emtr_trees.tree, 200)              AS tree_preview,
        CHAR_LENGTH(emtr_trees.tree)            AS tree_len,
        rep, iter, epsilon, n_leaves, n_edges, root, is_current,
        created_at, updated_at
      FROM emtr_trees
      ${where}
      ORDER BY created_at DESC
      LIMIT ? OFFSET ?
      `,
      [fullInt, fullInt, fullInt, job_id, limit, offset]
    );

    return NextResponse.json({ items: rows, trees: rows }, { status: 200 });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ error: msg }, { status: 500 });
  }
}

/* ---------------- POST ---------------- */

export async function POST(req: Request) {
  try {
    const ct = (req.headers.get("content-type") || "").toLowerCase();
    if (!ct.includes("application/json")) {
      return NextResponse.json({ error: "Content-Type must be application/json" }, { status: 415 });
    }

    const body = (await req.json().catch(() => null)) as Record<string, unknown> | null;
    if (!body) {
      return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
    }

    const job_id = String(body.job_id ?? "").trim();
    // Allow alias "newick" too
    const tree = String(body.tree ?? body.newick ?? "");

    if (!job_id || !tree) {
      return NextResponse.json({ error: "job_id and tree are required" }, { status: 400 });
    }

    const source = (body.source ?? "familyjoining") as string;
    const method = (body.method ?? "5foldcv") as string;
    const label = (body.label ?? "final") as string | null;

    const fmtRaw = String(body.format ?? "newick").toLowerCase();
    const format = (["newick", "nexus", "json"].includes(fmtRaw) ? fmtRaw : "newick") as
      | "newick" | "nexus" | "json";

    // Optional: tree with internal node labels (h_*)
    const tree_labeled = body.tree_labeled != null ? String(body.tree_labeled) : tree;

    // Edge list: accept array or JSON string
    const edge_list = normalizeEdgeList(body.edge_list);
    if (!edge_list) {
      return NextResponse.json({ error: "edge_list must be a JSON array of [u, v, length]" }, { status: 400 });
    }
    const edgeJson = JSON.stringify(edge_list);

    // Numerics
    const rep = toNullableNumber(body.rep);
    const iter = toNullableNumber(body.iter);
    const epsilon = toNullableNumber(body.epsilon);

    // Infer counts if not provided
    let n_leaves = toNullableNumber(body.n_leaves);
    const n_edges = toNullableNumber(body.n_edges) ?? edge_list.length;
    if (n_leaves == null) {
      n_leaves = inferLeafCountFromEdges(edge_list);
    }

    const root = body.root == null ? null : String(body.root);
    const is_current = coerceYesNo(body.is_current);

    // SHA-256 of the (stripped) tree text for dedupe/integrity
    const hashHex = crypto.createHash("sha256").update(tree).digest("hex");

    await q(
      `
      INSERT INTO emtr_trees
        (job_id, source, method, label, format,
         tree, tree_labeled, edge_list,
         rep, iter, epsilon, n_leaves, n_edges, root, is_current, tree_sha256)
      VALUES
        (?, ?, ?, ?, ?,
         ?, ?, CAST(? AS JSON),
         ?, ?, ?, ?, ?, ?, ?, UNHEX(?))
      `,
      [
        job_id, source, method, label, format,
        tree, tree_labeled, edgeJson,
        rep, iter, epsilon, n_leaves, n_edges, root, is_current, hashHex,
      ]
    );

    type IdRow = { [k: string]: unknown; id: number };
    const rows2 = await q<IdRow>(
      `
      SELECT id
      FROM emtr_trees
      WHERE job_id = ? AND tree_sha256 = UNHEX(?)
      ORDER BY id DESC
      LIMIT 1
      `,
      [job_id, hashHex]
    );
    const newId = rows2[0]?.id;

    if (!newId) {
      return NextResponse.json({ error: "inserted row not found" }, { status: 500 });
    }

    if (is_current) {
      await q(`UPDATE emtr_trees SET is_current=0 WHERE job_id=? AND id<>?`, [job_id, newId]);
    }

    return NextResponse.json({ id: newId, sha256: hashHex }, { status: 201 });
  } catch (e) {
    const msg = e instanceof Error ? e.message : "Unknown error";
    return NextResponse.json({ error: msg }, { status: 500 });
  }
}
