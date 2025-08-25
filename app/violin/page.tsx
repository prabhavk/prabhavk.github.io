// app/violin/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { buildExportPrefix } from "@/lib/exportName";
import dynamic from "next/dynamic";
import type { Data, Layout, PlotlyHTMLElement } from "plotly.js";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

type MethodName = "Parsimony" | "Dirichlet" | "SSH";

// Enforce this order on the x-axis and when creating traces
const METHOD_ORDER: MethodName[] = ["Dirichlet", "Parsimony", "SSH"];

// Nodes (dropdown)
const NODES = Array.from({ length: 17 }, (_, i) => `h_${21 + i}`);

type SeriesResp = {
  job_id: string;
  root?: string;
  series: { Parsimony: number[]; Dirichlet: number[]; SSH: number[] };
  counts: { Parsimony: number; Dirichlet: number; SSH: number };
};
type ApiErr = { error: string };

// ------- tiny guards -------
function hasProp<K extends string>(obj: unknown, key: K): obj is Record<K, unknown> {
  return typeof obj === "object" && obj !== null && key in obj;
}
function isErr(x: unknown): x is ApiErr {
  return hasProp(x, "error") && typeof x.error === "string";
}

export default function ViolinPage() {
  const [job, setJob] = useState<string>("");
  const [root, setRoot] = useState<string>("h_32");
  const [data, setData] = useState<SeriesResp | null>(null);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // ---- EXPORT / DOWNLOAD STATE & HELPERS ----
  const exportPrefix = useMemo(() => buildExportPrefix?.({ job }) ?? (job || "violin"), [job]);

  // all chart refs live here; key = your chosen id
  const violinRefs = useRef<Record<string, PlotlyHTMLElement | null>>({});

  // register a chart under a stable id
  const registerViolin = useCallback(
    (id: string) =>
      (_figure: unknown, gd: PlotlyHTMLElement) => {
        violinRefs.current[id] = gd;
      },
    []
  );

  // Plotly toImage access
  const toImage = useCallback(
    async (gd: PlotlyHTMLElement, opts: { format?: "png" | "svg" | "jpeg" | "webp"; scale?: number }) => {
      const w = window as unknown as { Plotly?: { toImage: (el: PlotlyHTMLElement, o: typeof opts) => Promise<string> } };
      if (!w.Plotly?.toImage) throw new Error("Plotly.toImage not available");
      return w.Plotly.toImage(gd, opts);
    },
    []
  );

  const downloadURI = (uri: string, name: string) => {
    const a = document.createElement("a");
    a.href = uri;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
  };

  // Compact fixed-width UTC stamp: YYYYMMDD_HHMMSSZ
  const utcStamp = () =>
    new Date().toISOString().replace(/\.\d{3}Z$/, "Z").replace(/[-:]/g, "").replace("T", "_");

  // Download a single chart by its id
  const downloadViolin = useCallback(
    async (id: string, fmt: "png" | "svg" = "png") => {
      const gd = violinRefs.current[id];
      if (!gd) return;
      const url = await toImage(gd, { format: fmt, scale: 2 });
      const fname = `${exportPrefix}__violin__${id}__t${utcStamp()}.${fmt}`;
      downloadURI(url, fname);
    },
    [toImage, exportPrefix]
  );

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

  // ----- helpers -----
  function median(arr: number[]): number {
    const a = arr.filter((v) => Number.isFinite(v)).slice().sort((x, y) => x - y);
    const n = a.length;
    if (!n) return NaN;
    return n % 2 ? a[(n - 1) >> 1] : (a[n / 2 - 1] + a[n / 2]) / 2;
  }

  // Build violin traces in required order — black on white
  const violinTraces: Partial<Data>[] = useMemo(() => {
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
      line: { color: "black" },
      fillcolor: "rgba(0,0,0,0.15)",
      opacity: 0.85,
      marker: { color: "black" },
    }));
  }, [data]);

  // Median line across the three categories (black)
  const medianLineTrace: Partial<Data> | null = useMemo(() => {
    if (!data) return null;
    const medians = METHOD_ORDER.map((m) => median(data.series[m]));
    return {
      type: "scatter",
      mode: "lines+markers",
      x: METHOD_ORDER as unknown as string[], // Plotly wants string[] here
      y: medians,
      line: { color: "black", width: 2 },
      marker: { color: "black", size: 8 },
      name: "Median",
      showlegend: false,
      hoverinfo: "skip" as const,
    };
  }, [data]);

  const layout: Partial<Layout> = {
    title: { text: `ll_final · root=${root}`, font: { color: "black" } },
    xaxis: {
      title: { text: "Method", font: { color: "black" } },
      tickfont: { color: "black" },
      categoryorder: "array",
      categoryarray: METHOD_ORDER,
      zerolinecolor: "#e5e7eb",
      gridcolor: "#f3f4f6",
    },
    yaxis: {
      title: { text: "ll_final", font: { color: "black" } },
      tickfont: { color: "black" },
      zerolinecolor: "#e5e7eb",
      gridcolor: "#f3f4f6",
    },
    margin: { l: 60, r: 20, t: 40, b: 50 },
    height: 480,
    showlegend: false,
    paper_bgcolor: "white",
    plot_bgcolor: "white",
  };

  return (
    <div className="p-6 max-w-[1200px] mx-auto bg-white text-black">
      <h1 className="text-2xl font-bold mb-4">Violin Plots of Maximum Log Likelihood</h1>

      <div className="flex flex-wrap gap-3 items-end mb-4">
        <label className="flex flex-col">
          <span className="text-sm">Root node</span>
          <select className="border rounded px-3 py-2" value={root} onChange={(e) => setRoot(e.target.value)}>
            {NODES.map((n) => (
              <option key={n} value={n}>
                {n}
              </option>
            ))}
          </select>
        </label>

        <div className="text-sm">
          {job ? (
            <>
              Job: <span className="font-mono">{job}</span>
            </>
          ) : (
            "No job selected"
          )}
        </div>

        <div className="ml-auto flex gap-2">
          {/* Downloads */}
          <button
            className="px-4 py-2 rounded bg-black text-white disabled:opacity-50"
            onClick={() => downloadViolin("ll_final_violin", "png")}
            disabled={!data}
            title="Download current violin as PNG"
          >
            Download PNG
          </button>          
        </div>
      </div>

      {loading ? (
        <div className="p-4 border rounded">Loading…</div>
      ) : err ? (
        <div className="p-4 border rounded text-red-600">{err}</div>
      ) : data ? (
        <>
          <div className="grid grid-cols-3 gap-3 mb-4">
            {METHOD_ORDER.map((m) => (
              <div key={m} className="p-3 border rounded bg-white flex items-center justify-between">
                <div className="font-semibold flex items-center gap-2">
                  <span className="inline-block w-3 h-3 rounded-full bg-black" aria-hidden />
                  <span>{m}</span>
                </div>
                <div className="text-sm">{data.counts[m]} repetitions</div>
              </div>
            ))}
          </div>

          {/* Violin chart + median line */}
          <div className="bg-white border rounded p-2">
            <Plot
              data={medianLineTrace ? [...violinTraces, medianLineTrace] : violinTraces}
              layout={layout}
              config={{ displayModeBar: false, responsive: true }}
              style={{ width: "100%", height: "100%" }}
              onInitialized={registerViolin("ll_final_violin")}
              onUpdate={registerViolin("ll_final_violin")}
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
