// app/violin/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import dynamic from "next/dynamic";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

type MethodName = "Parsimony" | "Dirichlet" | "SSH";

// Enforce this order everywhere
const METHOD_ORDER: MethodName[] = ["Dirichlet", "Parsimony", "SSH"];

// Nodes
const NODES = Array.from({ length: 17 }, (_, i) => `h_${21 + i}`);

type SeriesResp = {
  job_id: string;
  root?: string;
  series: { Parsimony: number[]; Dirichlet: number[]; SSH: number[] };
  counts: { Parsimony: number; Dirichlet: number; SSH: number };
};
type ApiErr = { error: string };

type ViolinTrace = {
  type: "violin";
  name: MethodName;
  y: number[];
  box?: { visible: boolean };
  meanline?: { visible: boolean };
  points?: "all" | "none" | "outliers" | "suspectedoutliers";
  jitter?: number;
  scalemode?: "width" | "count" | "area";
  showlegend?: boolean;
  line?: { color?: string };
  fillcolor?: string;
  opacity?: number;
  marker?: { color?: string };
};

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

  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId") ?? "";
      setJob(saved);
    } catch {/* ignore */}
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
        throw new Error(`Expected JSON, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`);
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

  useEffect(() => { if (job) void load(); }, [job, root, load]);

  const violinTraces: ViolinTrace[] = useMemo(() => {
    if (!data) return [];
    return METHOD_ORDER.map((m) => ({
      type: "violin",
      name: m,
      y: data.series[m],
      box: { visible: true },
      meanline: { visible: true },
      points: "all",
      jitter: 0.3,
      scalemode: "width",
      showlegend: false,
      line: { color: "#000" },               // black outlines
      fillcolor: "rgba(0,0,0,0.5)",          // black fill (semi-transparent)
      opacity: 0.8,
      marker: { color: "#000" },             // black points
    }));
  }, [data]);

  const medianLineTrace = useMemo(() => {
    if (!data) return null;
    const y = METHOD_ORDER.map((m) => median(data.series[m]));
    return {
      type: "scatter" as const,
      mode: "lines+markers" as const,
      x: METHOD_ORDER,
      y,
      line: { color: "black", width: 2 },
      marker: { color: "black", size: 6 },
      name: "median",
      showlegend: false,
      hoverinfo: "x+y",
    };
  }, [data]);

  return (
    <div className="p-6 max-w-[1200px] mx-auto">
      <h1 className="text-2xl font-bold mb-4">
        Violin Plots of <code>ll_final</code>
      </h1>

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
          {/* Method badges in desired order; all black */}
          <div className="grid grid-cols-3 gap-3 mb-4">
            {METHOD_ORDER.map((m) => (
              <div key={m} className="p-3 border rounded bg-white flex items-center justify-between">
                <div className="font-semibold flex items-center gap-2">
                  <span className="inline-block w-3 h-3 rounded-full" style={{ backgroundColor: "#000" }} aria-hidden />
                  <span style={{ color: "#000" }}>{m}</span>
                </div>
                <div className="text-sm" style={{ color: "#000" }}>
                  {data.counts[m]} repetitions
                </div>
              </div>
            ))}
          </div>

          {/* Violin chart with medians joined by black line */}
          <div className="bg-white border rounded p-2">
            <Plot
              data={medianLineTrace ? [...violinTraces, medianLineTrace] : violinTraces}
              layout={{
                title: { text: `ll_final · root=${root}` },
                xaxis: {
                  title: { text: "Method" },
                  categoryorder: "array",
                  categoryarray: METHOD_ORDER, // Dirichlet, Parsimony, SSH
                },
                yaxis: { title: { text: "ll_final" } },
                margin: { l: 60, r: 20, t: 40, b: 50 },
                height: 480,
                showlegend: false,
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

/** median of finite numbers; returns NaN if empty after filtering */
function median(arr: number[]): number {
  const a = arr.filter((v) => Number.isFinite(v)).slice().sort((x, y) => x - y);
  const n = a.length;
  if (!n) return NaN;
  const mid = Math.floor(n / 2);
  return n % 2 ? a[mid] : (a[mid - 1] + a[mid]) / 2;
}
