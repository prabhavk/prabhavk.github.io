// app/api/cake/route.ts
import { NextResponse } from "next/server";
import { query as q, type SQLParam } from "@/lib/db";
import crypto from "node:crypto";

export const runtime = "nodejs";
// If you want to always hit the DB and never cache at the edge:
// export const dynamic = "force-dynamic";

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
  tree?: string | null;
  tree_labeled?: string | null;
  edge_list?: Edge[] | string | null;
  tree_preview: string;
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

/* ---------- Optional columns discovery (cached) ---------- */

type SchemaInfo = { hasTreeLabeled: boolean; edgeCol: "edge_list" | "edge_list_json" | null };
let schemaPromise: Promise<SchemaInfo> | null = null;

async function detectSchema(): Promise<SchemaInfo> {
  if (schemaPromise) return schemaPromise;
  schemaPromise = (async () => {
    try {
      const rows = await q<{ column_name: string }>(
        `
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = DATABASE()
          AND table_name = 'emtr_trees'
          AND column_name IN ('tree_labeled','edge_list','edge_list_json')
        `
      );
      const names = new Set((rows ?? []).map(r => String(r.column_name)));
      const hasTreeLabeled = names.has("tree_labeled");
      const edgeCol = names.has("edge_list")
        ? "edge_list"
        : (names.has("edge_list_json") ? "edge_list_json" : null);
      return { hasTreeLabeled, edgeCol };
    } catch {
      return { hasTreeLabeled: false, edgeCol: null };
    }
  })();
  return schemaPromise;
}

/* ---------- Edge helpers ---------- */

type EdgeInput = Edge | [unknown, unknown, unknown];

function normalizeEdgeList(raw: unknown): Edge[] | null {
  if (raw == null) return null;
  let x: unknown = raw;

  if (typeof x === "string") {
    const t = x.trim();
    if (!t) return null;
    try {
      x = JSON.parse(t);
    } catch {
      return null;
    }
  }

  if (!Array.isArray(x)) return null;

  const out: Edge[] = [];
  for (const e of x as EdgeInput[]) {
    if (!Array.isArray(e) || e.length !== 3) return null;
    const a = e[0], b = e[1], w = Number(e[2]);
    if (typeof a !== "string" || typeof b !== "string" || !Number.isFinite(w)) return null;
    out.push([a, b, w]);
  }
  return out;
}

function inferLeafCountFromEdges(edges: Edge[]): number {
  const names = new Set<string>();
  for (const [u, v] of edges) {
    names.add(u);
    names.add(v);
  }
  let cnt = 0;
  for (const name of names) if (!/^h_/.test(name)) cnt++;
  return cnt;
}

/* ---------------- GET ---------------- */

export async function GET(req: Request) {
  try {
    const url = new URL(req.url);

    const job_id = (url.searchParams.get("job_id") || url.searchParams.get("job") || "").trim();
    const currentOnly = truthy(url.searchParams.get("current"));
    const full = truthy(url.searchParams.get("full"));

    const limit = clampInt(url.searchParams.get("limit"), 50, 1, 200);
    const offset = clampInt(url.searchParams.get("offset"), 0, 0, 1_000_000);

    if (!job_id) {
      return NextResponse.json({ error: "job_id required" }, { status: 400 });
    }

    const { hasTreeLabeled, edgeCol } = await detectSchema();

    const cols: string[] = [
      "id", "job_id", "source", "method", "label", "format",
      full ? "emtr_trees.tree AS tree" : "NULL AS tree",
      "LEFT(emtr_trees.tree, 200) AS tree_preview",
      "CHAR_LENGTH(emtr_trees.tree) AS tree_len",
      "rep", "iter", "epsilon", "n_leaves", "n_edges", "root", "is_current",
      "created_at", "updated_at",
    ];

    // Insert optional columns in stable positions
    cols.splice(7, 0, hasTreeLabeled
      ? (full ? "emtr_trees.tree_labeled AS tree_labeled" : "NULL AS tree_labeled")
      : "NULL AS tree_labeled");

    cols.splice(8, 0, edgeCol
      ? (full ? `emtr_trees.${edgeCol} AS edge_list` : "NULL AS edge_list")
      : "NULL AS edge_list");

    const where = `WHERE job_id = ? ${currentOnly ? "AND is_current=1" : ""}`;
    const sql = `
      SELECT ${cols.join(", ")}
      FROM emtr_trees
      ${where}
      ORDER BY created_at DESC
      LIMIT ? OFFSET ?
    `;

    const rows = await q<TreeRow>(sql, [job_id, limit, offset]);

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
    const tree = String(body.tree ?? body.newick ?? "");

    if (!job_id || !tree) {
      return NextResponse.json({ error: "job_id and tree are required" }, { status: 400 });
    }

    const source = String(body.source ?? "familyjoining");
    const method = String(body.method ?? "5foldcv");
    const label = (body.label == null ? null : String(body.label)) as string | null;

    const fmtRaw = String(body.format ?? "newick").toLowerCase();
    const format = (["newick", "nexus", "json"].includes(fmtRaw) ? fmtRaw : "newick") as
      | "newick" | "nexus" | "json";

    const { hasTreeLabeled, edgeCol } = await detectSchema();

    const tree_labeled = body.tree_labeled != null ? String(body.tree_labeled) : null;

    const edge_list = normalizeEdgeList(body.edge_list);
    const edgeJson = edge_list ? JSON.stringify(edge_list) : null;

    const rep = toNullableNumber(body.rep);
    const iter = toNullableNumber(body.iter);
    const epsilon = toNullableNumber(body.epsilon);

    let n_leaves = toNullableNumber(body.n_leaves);
    let n_edges = toNullableNumber(body.n_edges);
    if (edge_list) {
      n_edges = n_edges ?? edge_list.length;
      n_leaves = n_leaves ?? inferLeafCountFromEdges(edge_list);
    }

    const root = body.root == null ? null : String(body.root);
    const is_current = coerceYesNo(body.is_current);

    const hashHex = crypto.createHash("sha256").update(tree).digest("hex");

    // ---- Build INSERT dynamically
    const cols: string[] = [
      "job_id", "source", "method", "label", "format", "tree",
    ];
    const placeholders: string[] = ["?","?","?","?","?","?"];
    const params: SQLParam[] = [
      job_id, source, method, (label ?? null), format, tree,
    ];

    if (hasTreeLabeled) {
      cols.push("tree_labeled");
      placeholders.push("?");
      params.push((tree_labeled ?? tree) as SQLParam);
    }

    if (edgeCol && edgeJson !== null) {
      cols.push(edgeCol);
      placeholders.push("?");
      params.push(edgeJson);
    }

    cols.push("rep","iter","epsilon","n_leaves","n_edges","root","is_current","tree_sha256");
    placeholders.push("?","?","?","?","?","?","?","UNHEX(?)");
    params.push(rep, iter, epsilon, n_leaves, n_edges, (root ?? null), is_current, hashHex);

    const sql = `
      INSERT INTO emtr_trees
        (${cols.join(", ")})
      VALUES
        (${placeholders.join(", ")})
    `;

    await q(sql, params);

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
