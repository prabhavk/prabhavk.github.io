// app/violin/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import dynamic from "next/dynamic";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

// Plotly default categorical palette (ensures labels match trace colors)
const PLOTLY_COLORWAY = [
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
] as const;

type MethodName = "Parsimony" | "Dirichlet" | "SSH";
const METHOD_INDEX: Record<MethodName, number> = { Parsimony: 0, Dirichlet: 1, SSH: 2 };
const methodColor = (m: MethodName) => PLOTLY_COLORWAY[METHOD_INDEX[m]];

// Nodes
const NODES = Array.from({ length: 17 }, (_, i) => `h_${21 + i}`);

type SeriesResp = {
  job_id: string;
  root?: string;
  series: { Parsimony: number[]; Dirichlet: number[]; SSH: number[] };
  counts: { Parsimony: number; Dirichlet: number; SSH: number };
};
type ApiErr = { error: string };

// minimal trace typing
type ViolinTrace = {
  type: "violin";
  name: MethodName;
  y: number[];
  box?: { visible: boolean };
  meanline?: { visible: boolean };
  points?: "all" | "none" | "outliers" | "suspectedoutliers";
  jitter?: number;
  scalemode?: "width" | "count" | "area";
};

// ------- tiny guards -------
function hasProp<K extends string>(obj: unknown, key: K): obj is Record<K, unknown> {
  return typeof obj === "object" && obj !== null && key in obj;
}
function isErr(x: unknown): x is ApiErr {
  return hasProp(x, "error") && typeof x.error === "string";
}

export default function ViolinPage() {
  const [job, setJob] = useState<string>("");
  const [root, setRoot] = useState<string>(NODES[0]);
  const [data, setData] = useState<SeriesResp | null>(null);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // Read the job chosen on the Precomputed Results page
  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId") ?? "";
      setJob(saved);
    } catch {
      /* ignore */
    }
  }, []);

  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Please pick a job on the Precomputed Results page.");
      setData(null);
      return;
    }
    setLoading(true);
    setErr(null);
    try {
      const u = new URL("/api/violin", window.location.origin);
      u.searchParams.set("job", job);
      if (root) u.searchParams.set("root", root);
      const res = await fetch(u.toString(), { cache: "no-store" });

      const ct = res.headers.get("content-type") || "";
      if (!ct.includes("application/json")) {
        const body = await res.text().catch(() => "");
        throw new Error(
          `Expected JSON from ${u.pathname}, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`
        );
      }

      const json: unknown = await res.json();
      if (!res.ok || isErr(json)) throw new Error(isErr(json) ? json.error : `HTTP ${res.status}`);
      setData(json as SeriesResp);
    } catch (e: unknown) {
      setErr(e instanceof Error ? e.message : "Failed to load");
      setData(null);
    } finally {
      setLoading(false);
    }
  }, [job, root]);

  // Auto-load when job/root change
  useEffect(() => {
    if (job) void load();
  }, [job, root, load]);

  const traces: ViolinTrace[] = useMemo(() => {
    if (!data) return [];
    // No explicit colors: Plotly will use its default colorway (matching labels we render)
    return (["Parsimony", "Dirichlet", "SSH"] as MethodName[]).map((m) => ({
      type: "violin",
      name: m,
      y: data.series[m],
      box: { visible: true },
      meanline: { visible: true },
      points: "all",
      jitter: 0.3,
      scalemode: "width",
    }));
  }, [data]);

  return (
    <div className="p-6 max-w-[1200px] mx-auto">
      <h1 className="text-2xl font-bold mb-4">Violin Plots of <code>ll_final</code></h1>

      <div className="flex flex-wrap gap-3 items-end mb-4">
        <label className="flex flex-col">
          <span className="text-sm">Root node</span>
          <select className="border rounded px-3 py-2" value={root} onChange={(e) => setRoot(e.target.value)}>
            {NODES.map((n) => (
              <option key={n} value={n}>{n}</option>
            ))}
          </select>
        </label>

        <div className="text-sm text-white">
          {job ? <>Job: <span className="font-mono text-white">{job}</span></> : "No job selected"}
        </div>

        <button
          className="ml-auto px-4 py-2 rounded bg-black text-white disabled:opacity-50"
          onClick={() => void load()}
          disabled={loading || !job}
        >
          {loading ? "Loading…" : "Reload"}
        </button>
      </div>

      {loading ? (
        <div className="p-4 border rounded">Loading…</div>
      ) : err ? (
        <div className="p-4 border rounded text-red-600">{err}</div>
      ) : data ? (
        <>
          {/* Quick counts with badges colored to Plotly's default hues */}
          <div className="grid grid-cols-3 gap-3 mb-4">
            {(["Parsimony","Dirichlet","SSH"] as const).map((m) => (
              <div key={m} className="p-3 border rounded bg-white flex items-center justify-between">
                <div className="font-semibold flex items-center gap-2">
                  <span
                    className="inline-block w-3 h-3 rounded-full"
                    style={{ backgroundColor: methodColor(m) }}
                    aria-hidden
                  />
                  <span style={{ color: methodColor(m) }}>{m}</span>
                </div>
                <div className="text-sm" style={{ color: methodColor(m) }}>
                  {data.counts[m]} repetitions
                </div>
              </div>
            ))}
          </div>

          {/* Violin chart */}
          <div className="bg-white border rounded p-2">
            <Plot
              data={traces}
              layout={{
                title: { text: `ll_final · root=${root}` },
                colorway: PLOTLY_COLORWAY as unknown as string[], // keep mapping stable
                xaxis: { title: { text: "Method" } },
                yaxis: { title: { text: "ll_final" } },
                margin: { l: 60, r: 20, t: 40, b: 50 },
                height: 480,
              }}
              config={{ displayModeBar: false, responsive: true }}
              style={{ width: "100%", height: "100%" }}
            />
          </div>
        </>
      ) : (
        <div className="p-4 border rounded">
          No data. {job ? "Try Reload." : "Select a job on the Precomputed Results page."}
        </div>
      )}
    </div>
  );
}
