"use client";

import TreeViewer from "@/components/tree-viewer";

type YesNo = 0 | 1;

export type RowLite = {
  id: number;
  job_id: string;
  label: string | null;
  method: string;
  source: string;
  format: string | null;
  rep: number | null;
  iter: number | null;
  epsilon: number | null;
  n_leaves: number | null;
  n_edges: number | null;
  root_name: string | null;
  is_current: YesNo | boolean;
  created_at: string;
  updated_at: string;
};

type Props = {
  row: RowLite;
  /** Newick string already normalized/semicolon-terminated by the server page */
  newick: string;
};

export default function TreeDetailClient({ row, newick }: Props) {
  const created = row.created_at ? new Date(String(row.created_at)).toLocaleString() : "";

  return (
    <div className="space-y-6">
      <div className="grid gap-2 text-sm text-gray-800">
        <div>
          <span className="font-medium">Method:</span> {row.method}{" "}
          <span className="ml-4 font-medium">Source:</span> {row.source}{" "}
          <span className="ml-4 font-medium">Format:</span> {row.format ?? "—"}
        </div>
        <div>
          <span className="font-medium">Current:</span>{" "}
          {row.is_current ? "Yes" : "No"}{" "}
          <span className="ml-4 font-medium">Created:</span> {created}
        </div>
        <div>
          <span className="font-medium">rep:</span> {row.rep ?? "—"}{" "}
          <span className="ml-4 font-medium">iter:</span> {row.iter ?? "—"}{" "}
          <span className="ml-4 font-medium">ε:</span> {row.epsilon ?? "—"}{" "}
          <span className="ml-4 font-medium">leaves:</span> {row.n_leaves ?? "—"}{" "}
          <span className="ml-4 font-medium">edges:</span> {row.n_edges ?? "—"}
        </div>
      </div>

      {newick ? (
        <>
          <div>
            <div className="text-sm font-medium mb-1">Viewer</div>
            <TreeViewer newick={newick} height={480} />
          </div>

          <div>
            <div className="text-sm font-medium mb-1">Newick</div>
            <pre className="whitespace-pre-wrap break-words border rounded p-3 bg-gray-50">
              {newick}
            </pre>
          </div>
        </>
      ) : (
        <div className="text-red-600 text-sm">No tree text stored for this row.</div>
      )}
    </div>
  );
}
