// app/api/trees/route.ts
import { NextResponse } from "next/server";
import { query as q } from "@/lib/db";
import crypto from "node:crypto";

export const runtime = "nodejs"; // ensure node runtime so `node:crypto` works

type YesNo = 0 | 1;

type TreeRow = {
  // Satisfy lib/db generic constraint
  [key: string]: unknown;

  id: number;
  job_id: string;
  source: string;
  method: string;
  label: string | null;
  format: "newick" | "nexus" | "json" | string;
  tree?: string | null;     // present when full=1 (or null otherwise)
  tree_preview: string;     // always present
  tree_len: number;
  rep: number | null;
  iter: number | null;
  epsilon: number | null;
  n_leaves: number | null;
  n_edges: number | null;
  root_name: string | null;
  is_current: YesNo;
  created_at: string;
  updated_at: string;
};

export async function GET(req: Request) {
  const url = new URL(req.url);

  // Accept both ?job_id= and ?job=
  const job_id = url.searchParams.get("job_id") || url.searchParams.get("job");
  const currentOnly = url.searchParams.get("current") === "1";
  const full = url.searchParams.get("full") === "1" || url.searchParams.get("full") === "true";
  const fullInt = full ? 1 : 0;

  const limit = Math.min(parseInt(url.searchParams.get("limit") || "50", 10), 200);
  const offset = Math.max(parseInt(url.searchParams.get("offset") || "0", 10), 0);

  if (!job_id) {
    return NextResponse.json({ error: "job_id required" }, { status: 400 });
  }

  // Always return both: `tree` (nullable unless full=1) and `tree_preview`
  // Using IF(?=1, tree, NULL) keeps the result shape consistent for clients.
  const where = `WHERE job_id = ? ${currentOnly ? "AND is_current=1" : ""}`;
  const rows = await q<TreeRow>(
    `SELECT id, job_id, source, method, label, format,
            IF(?=1, tree, NULL) AS tree,
            LEFT(tree, 200) AS tree_preview,
            CHAR_LENGTH(tree) AS tree_len,
            rep, iter, epsilon, n_leaves, n_edges, root_name, is_current,
            created_at, updated_at
     FROM emtr_trees
     ${where}
     ORDER BY created_at DESC
     LIMIT ? OFFSET ?`,
    [fullInt, job_id, limit, offset]
  );

  return NextResponse.json({ items: rows, trees: rows }, { status: 200 });
}

export async function POST(req: Request) {
  const body = await req.json().catch(() => null);
  if (!body?.job_id || !body?.tree) {
    return NextResponse.json({ error: "job_id and tree are required" }, { status: 400 });
  }

  const {
    job_id,
    tree,
    source = "familyjoining",
    method = "5foldcv",   // keep your default as requested
    label = "final",
    format = "newick",
    rep = null,
    iter = null,
    epsilon = null,
    n_leaves = null,
    n_edges = null,
    root_name = null,
    is_current = 0 as YesNo,
  } = body as {
    job_id: string;
    tree: string;
    source?: string;
    method?: string;
    label?: string | null;
    format?: string;
    rep?: number | null;
    iter?: number | null;
    epsilon?: number | null;
    n_leaves?: number | null;
    n_edges?: number | null;
    root_name?: string | null;
    is_current?: YesNo;
  };

  // Normalize + hash (dedupe/integrity). Keep exactly as stored.
  const hashHex = crypto.createHash("sha256").update(String(tree)).digest("hex");

  // Insert
  await q(
    `INSERT INTO emtr_trees
     (job_id, source, method, label, format, tree, rep, iter, epsilon, n_leaves, n_edges, root_name, is_current, tree_sha256)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, UNHEX(?))`,
    [job_id, source, method, label, format, tree, rep, iter, epsilon, n_leaves, n_edges, root_name, is_current, hashHex]
  );

  // Fetch the inserted id robustly without relying on session-bound LAST_INSERT_ID()
  // We use the content-hash + job_id to find the most recent insert.
  type IdRow = { [key: string]: unknown; id: number };
  const rows2 = await q<IdRow>(
    `SELECT id FROM emtr_trees
     WHERE job_id = ? AND tree_sha256 = UNHEX(?)
     ORDER BY id DESC
     LIMIT 1`,
    [job_id, hashHex]
  );
  const newId = rows2[0]?.id;
  if (!newId) {
    return NextResponse.json({ error: "inserted row not found" }, { status: 500 });
  }

  // If this one is current, clear other currents for the same job
  if (is_current) {
    await q(`UPDATE emtr_trees SET is_current=0 WHERE job_id=? AND id<>?`, [job_id, newId]);
  }

  return NextResponse.json({ id: newId }, { status: 201 });
}
