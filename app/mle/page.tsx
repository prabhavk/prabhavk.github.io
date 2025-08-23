// app/mle/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import { ProbHists } from "@/components/ProbHists";

/** API payload */
type MleOk = {
  job_id?: string;
  method?: string;
  root_name?: string | null;   // used to exclude the root matrix in ProbHists
  root: number[] | null;       // expected length 4
  trans: unknown;              // 4x4 array OR Record<string, 4x4 array> OR flat 16
  D_pi?: number[] | null;
  D_M?: number[] | null;       // Dirichlet for transition rows (length 4)
};
type MleResp = MleOk | { error: string };

const SELECTED_JOB_KEY = "emtr:selectedJobId";

/** ---------- helpers ---------- */
function is4x4(x: unknown): x is number[][] {
  return (
    Array.isArray(x) &&
    x.length === 4 &&
    x.every(
      (r) =>
        Array.isArray(r) &&
        r.length === 4 &&
        r.every((v) => typeof v === "number" && Number.isFinite(v))
    )
  );
}

function flat16to4x4(flat: number[]): number[][] {
  return [
    flat.slice(0, 4),
    flat.slice(4, 8),
    flat.slice(8, 12),
    flat.slice(12, 16),
  ];
}

/** Convert various shapes into a 4x4 number matrix if possible. */
function to4x4(x: unknown): number[][] | null {
  if (is4x4(x)) return x;
  if (Array.isArray(x) && x.length === 16 && x.every((v) => Number.isFinite(v as number))) {
    return flat16to4x4(x as number[]);
  }
  return null;
}

/** Normalize API 'trans' into a map of key->4x4 matrices */
function normalizeTransToMap(trans: unknown): Record<string, number[][]> {
  const out: Record<string, number[][]> = {};

  // Single matrix (4x4 or flat16)
  const single = to4x4(trans);
  if (single) {
    out["1"] = single;
    return out;
  }

  // Object of many matrices (keys "1","2",...,"h_21", ...)
  if (trans && typeof trans === "object" && trans !== null) {
    for (const [k, v] of Object.entries(trans as Record<string, unknown>)) {
      const m = to4x4(v);
      if (m) out[k] = m;
    }
  }

  return out;
}

/** ---------- page ---------- */
export default function MLEinDPage() {
  const [job, setJob] = useState<string>("");
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  const [root, setRoot] = useState<number[] | null>(null);
  const [transMap, setTransMap] = useState<Record<string, number[][]> | null>(null);
  const [rootKey, setRootKey] = useState<string | null>(null);
  const [alphaPi, setAlphaPi] = useState<number[] | null>(null);
  const [alphaM, setAlphaM] = useState<number[] | null>(null); // ✅ NEW

  // pick job id from localStorage
  useEffect(() => {
    try {
      const saved = localStorage.getItem(SELECTED_JOB_KEY);
      if (saved) setJob(saved);
    } catch {
      /* ignore */
    }
  }, []);

  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Choose a job on the Precomputed Results page.");
      setRoot(null);
      setTransMap(null);
      setRootKey(null);
      setAlphaPi(null);
      setAlphaM(null);
      return;
    }

    setLoading(true);
    setErr(null);
    setRoot(null);
    setTransMap(null);
    setRootKey(null);
    setAlphaPi(null);
    setAlphaM(null);

    try {
      const u = new URL("/api/mle", window.location.origin);
      u.searchParams.set("job_id", job);
      u.searchParams.set("method", "dirichlet");

      const res = await fetch(u.toString(), { cache: "no-store" });
      const ct = res.headers.get("content-type") || "";
      if (!ct.toLowerCase().includes("application/json")) {
        const body = await res.text().catch(() => "");
        throw new Error(`Expected JSON, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`);
      }

      const j = (await res.json()) as MleResp;
      if (!res.ok || "error" in j) {
        throw new Error(("error" in j && j.error) || `HTTP ${res.status}`);
      }
      const data = j as MleOk;

      // Root probs
      const rootArr = Array.isArray(data.root) ? data.root.map(Number) : null;
      setRoot(rootArr);
      if (!rootArr || rootArr.length !== 4) {
        setErr((prev) => prev ?? "No root_prob_final (expected 4 values).");
      }

      // Transition matrices -> map
      const map = normalizeTransToMap(data.trans);
      setTransMap(Object.keys(map).length ? map : null);

      // root_name is optional; use it if present so ProbHists can exclude the root matrix
      const rKey = data.root_name && typeof data.root_name === "string" ? data.root_name : null;
      setRootKey(rKey);

      if (!Object.keys(map).length) {
        setErr((prev) => prev ?? "No transition probabilities found in trans_prob_final (expected 4×4).");
      }

      // Dirichlet alphas
      const aPi = Array.isArray(data.D_pi) && data.D_pi.length === 4 ? data.D_pi.map(Number) : null;
      setAlphaPi(aPi);

      const aM  = Array.isArray(data.D_M)  && data.D_M.length  === 4 ? data.D_M.map(Number)  : null;
      setAlphaM(aM); // ✅ NEW

    } catch (e) {
      setErr(e instanceof Error ? e.message : "Failed to load");
    } finally {
      setLoading(false);
    }
  }, [job]);

  useEffect(() => {
    if (job) void load();
  }, [job, load]);

  const rootOK = useMemo(() => Array.isArray(root) && root.length === 4, [root]);
  const transOK = useMemo(() => !!transMap && Object.keys(transMap).length > 0, [transMap]);

  return (
    <div className="p-6 max-w-[1400px] mx-auto space-y-4">
      <div className="flex flex-wrap items-end gap-3">
        <h1 className="text-2xl font-bold">Maximum likelihood estimate of GMM and density plot of Dirichlet parameters</h1>
        <div className="text-sm ml-4">
          Job: <span className="font-mono">{job || "(none)"} </span>
        </div>
        <button
          type="button"
          onClick={() => void load()}
          className="ml-auto px-4 py-2 rounded bg-black text-white disabled:opacity-50"
          disabled={!job || loading}
          title="Reload"
        >
          {loading ? "Loading…" : "Reload"}
        </button>
      </div>

      {loading ? (
        <div className="p-4 border rounded">Loading…</div>
      ) : err ? (
        <div className="p-4 border rounded text-red-600">{err}</div>
      ) : !(rootOK || transOK) ? (
        <div className="p-4 border rounded">No Dirichlet best repetition data available for this job.</div>
      ) : (
        <div className="space-y-4">
          <ProbHists
            root={rootOK ? root : null}
            // Pass the ENTIRE map so ProbHists can aggregate across matrices
            trans={transOK ? (transMap as unknown as Record<string, unknown>) : null}
            // Let ProbHists exclude the root matrix when aggregating
            rootKey={rootKey ?? undefined}
            rootName={rootKey ?? undefined}
            alphaPi={alphaPi}
            alphaM={alphaM}
          />
        </div>
      )}
    </div>
  );
}
