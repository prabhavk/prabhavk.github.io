// app/trees/[id]/page.tsx
import Link from "next/link";
import { query as q } from "@/lib/db";
import TreeDetailClient from "./TreeDetailClient"; // <- client component

export const dynamic = "force-dynamic";

type YesNo = 0 | 1;

interface TreeRow {
  [key: string]: unknown; // satisfy lib/db generic
  id: number;
  job_id: string;
  label: string | null;
  method: string;
  source: string;
  format: string | null;
  tree: string | null;
  rep: number | null;
  iter: number | null;
  epsilon: number | null;
  n_leaves: number | null;
  n_edges: number | null;
  root_name: string | null;
  is_current: YesNo;
  created_at: string | Date;
  updated_at: string | Date;
}

function ensureString(x: unknown): string {
  return typeof x === "string" ? x : x == null ? "" : String(x);
}
function ensureSemicolon(s: string): string {
  const t = s.trim();
  return t && !t.endsWith(";") ? `${t};` : t;
}

export default async function Page({
  // NOTE: params is a Promise here to satisfy your project's PageProps constraint
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const { id } = await params;
  const idNum = Number(id);

  if (!Number.isFinite(idNum)) {
    return (
      <div className="mx-auto max-w-5xl p-6">
        <p className="text-red-600">Invalid id.</p>
        <div className="mt-4">
          <Link href="/trees" className="text-blue-600 underline">
            ← Back to Trees
          </Link>
        </div>
      </div>
    );
  }

  const rows = await q<TreeRow>(
    `SELECT id, job_id, label, method, source, format, tree,
            rep, iter, epsilon, n_leaves, n_edges, root_name,
            is_current, created_at, updated_at
     FROM emtr_trees
     WHERE id = ?
     LIMIT 1`,
    [idNum]
  );

  const row = rows[0];
  if (!row) {
    return (
      <div className="mx-auto max-w-5xl p-6">
        <p className="text-gray-600">Tree not found.</p>
        <div className="mt-4">
          <Link href="/trees" className="text-blue-600 underline">
            ← Back to Trees
          </Link>
        </div>
      </div>
    );
  }

  // Normalize Newick text for the viewer
  const raw = ensureString(row.tree);
  const isNewick = String(row.format || "").toLowerCase() === "newick";
  const newick = isNewick ? ensureSemicolon(raw) : raw;

  return (
    <TreeDetailClient
      row={{
        id: row.id,
        job_id: ensureString(row.job_id),
        label: row.label,
        method: row.method,
        source: row.source,
        format: row.format,
        rep: row.rep,
        iter: row.iter,
        epsilon: row.epsilon,
        n_leaves: row.n_leaves,
        n_edges: row.n_edges,
        root_name: row.root_name,
        is_current: row.is_current,
        created_at: ensureString(row.created_at),
        updated_at: ensureString(row.updated_at),
      }}
      newick={newick}
    />
  );
}
