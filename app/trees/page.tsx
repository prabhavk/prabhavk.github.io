// app/trees/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import TreeViewer from "@/components/tree-viewer";

/* -------------------- Types -------------------- */
type YesNo = 0 | 1;

export type TreeRec = {
  id: number;
  job_id: string;
  source: string;
  method: string;
  label: string | null;
  format: "newick" | string;
  tree?: string | null;      // present when full=1
  newick?: string | null;    // some APIs use this field name
  is_current: YesNo | boolean;
  created_at?: string;
};

type ApiTrees = {
  ok?: boolean;
  error?: unknown;
  trees?: unknown;
  items?: unknown;
};

/* -------------------- Utilities -------------------- */
const isRecord = (v: unknown): v is Record<string, unknown> =>
  v !== null && typeof v === "object";

const isYesNoBool = (v: unknown): v is YesNo | boolean =>
  typeof v === "boolean" || v === 0 || v === 1;

const isOptString = (v: unknown): v is string | null | undefined =>
  v == null || typeof v === "string";

const isTreeRec = (v: unknown): v is TreeRec => {
  if (!isRecord(v)) return false;
  const id = v["id"];
  const job_id = v["job_id"];
  const source = v["source"];
  const method = v["method"];
  const format = v["format"];
  const is_current = v["is_current"];
  const label = v["label"];
  const tree = v["tree"];
  const newick = v["newick"];
  const created_at = v["created_at"];

  return (
    typeof id === "number" &&
    typeof job_id === "string" &&
    typeof source === "string" &&
    typeof method === "string" &&
    typeof format === "string" &&
    isYesNoBool(is_current) &&
    isOptString(label) &&
    isOptString(tree) &&
    isOptString(newick) &&
    isOptString(created_at)
  );
};

function getQueryParam(name: string): string {
  if (typeof window === "undefined") return "";
  const u = new URL(window.location.href);
  return u.searchParams.get(name) ?? "";
}

function normalizeNewick(n: unknown): string {
  if (typeof n !== "string") return "";
  const t = n.trim();
  return t && !t.endsWith(";") ? `${t};` : t;
}

/** Sum of all branch lengths in a Newick string (unrooted/rooted agnostic). */
function sumBranchLengths(nw: string): number {
  let total = 0;
  // numbers immediately after a colon, including scientific notation
  const re = /:\s*([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?)/g;
  let m: RegExpExecArray | null;
  while ((m = re.exec(nw)) !== null) {
    const v = Number(m[1]);
    if (Number.isFinite(v)) total += v;
  }
  return total;
}

/* -------------------- Data fetch -------------------- */
async function fetchCurrentTree(jobId: string, signal?: AbortSignal): Promise<TreeRec | null> {
  const u = new URL("/api/trees", window.location.origin);
  u.searchParams.set("job", jobId);
  u.searchParams.set("current", "1");
  u.searchParams.set("full", "1");

  const res = await fetch(u.toString(), { cache: "no-store", signal });
  const ct = (res.headers.get("content-type") || "").toLowerCase();
  const isJson = ct.includes("application/json");

  // Read body once
  const bodyText = await res.text();
  const bodyJson: unknown = isJson
    ? (() => { try { return JSON.parse(bodyText); } catch { return {}; } })()
    : {};

  if (!res.ok) {
    let msg = `HTTP ${res.status}`;
    if (isJson && isRecord(bodyJson)) {
      const errVal = (bodyJson as ApiTrees).error;
      if (typeof errVal === "string" && errVal.trim()) msg = errVal;
    } else if (bodyText) {
      msg = bodyText.slice(0, 200);
    }
    throw new Error(msg);
  }

  const payload: ApiTrees = isRecord(bodyJson) ? (bodyJson as ApiTrees) : {};
  const rowsUnknown: unknown[] = Array.isArray(payload.trees)
    ? (payload.trees as unknown[])
    : Array.isArray(payload.items)
    ? (payload.items as unknown[])
    : [];

  const rows: TreeRec[] = rowsUnknown.filter(isTreeRec);
  return rows[0] ?? null;
}

/* -------------------- Page -------------------- */
export default function TreesPage() {
  const [job, setJob] = useState<string>("");
  const [tree, setTree] = useState<TreeRec | null>(null);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // Derive job from URL or localStorage (no UI controls on this page)
  useEffect(() => {
    try {
      const fromUrl = getQueryParam("job") || getQueryParam("job_id");
      const fromLs = localStorage.getItem("emtr:selectedJobId") || "";
      setJob(fromUrl || fromLs || "");
    } catch {
      /* ignore */
    }
  }, []);

  // React to selection changes (e.g., other pages updating localStorage)
  useEffect(() => {
    function onStorage(ev: StorageEvent) {
      if (ev.key === "emtr:selectedJobId") {
        setJob(ev.newValue || "");
      }
    }
    window.addEventListener("storage", onStorage);
    return () => window.removeEventListener("storage", onStorage);
  }, []);

  const reload = useCallback(async () => {
    if (!job) {
      setTree(null);
      setErr("No job selected. Choose a job on the Precomputed Results page.");
      return;
    }
    setLoading(true);
    setErr(null);

    const ctrl = new AbortController();
    try {
      const t = await fetchCurrentTree(job, ctrl.signal);
      setTree(t);
      if (!t) setErr(`No current tree found for job ${job}.`);
    } catch (e) {
      setErr(e instanceof Error ? e.message : "Failed to load tree");
      setTree(null);
    } finally {
      setLoading(false);
    }
    // No need to abort here; provided as pattern for future concurrent reloads
    // return () => ctrl.abort();
  }, [job]);

  useEffect(() => {
    if (job) void reload();
  }, [job, reload]);

  // prefer `tree`, fall back to `newick`
  const newickRaw = tree?.tree ?? tree?.newick ?? "";
  const newick = useMemo(() => normalizeNewick(newickRaw), [newickRaw]);

  // Compute total tree length (sum of branch lengths)
  const totalLen = useMemo(() => {
    if (!newick) return null;
    const s = sumBranchLengths(newick);
    return Number.isFinite(s) ? s : null;
  }, [newick]);

  return (
    <div className="mx-auto max-w-6xl p-6 space-y-4">
      <div className="flex items-center gap-3">
        <h1 className="text-2xl font-semibold">Trees</h1>
        <div className="text-sm text-gray-600">
          Job: <span className="font-mono">{job || "(none)"}</span>
        </div>
        <button
          type="button"
          onClick={() => void reload()}
          className="ml-auto px-3 py-1.5 rounded bg-black text-white disabled:opacity-50 text-sm"
          disabled={!job || loading}
        >
          {loading ? "Loading…" : "Refresh"}
        </button>
      </div>

      {loading ? (
        <div className="p-3 border rounded">Loading…</div>
      ) : err ? (
        <div className="p-3 border border-red-600 rounded text-red-600">{err}</div>
      ) : !tree ? (
        <div className="p-3 border rounded">No current tree to display.</div>
      ) : (
        <>
          <div className="text-sm text-gray-600">
            <b>{tree.label ?? "Current tree"}</b> · method: <b>{tree.method}</b> · source:{" "}
            <b>{tree.source}</b>{" "}
            {tree.created_at ? <>· created: {new Date(tree.created_at).toLocaleString()}</> : null}
            {totalLen !== null ? (
              <> · total length (sum of branch lengths): <b>{totalLen.toFixed(4)}</b></>
            ) : null}
          </div>
          <TreeViewer newick={newick} height={520} />
        </>
      )}
    </div>
  );
}
