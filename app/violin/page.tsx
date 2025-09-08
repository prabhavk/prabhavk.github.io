// app/violin/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { buildExportPrefix } from "@/lib/exportName";
import dynamic from "next/dynamic";
import type { Data, Layout, PlotlyHTMLElement } from "plotly.js";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

type MethodName = "Parsimony" | "Dirichlet" | "HSS";
const METHOD_ORDER: MethodName[] = ["Dirichlet", "Parsimony", "HSS"];

// ----- Colors -----
const COLORS: Record<MethodName, string> = {
  Dirichlet: "#D2691E", // chocolate
  Parsimony: "#FF6B3D",
  HSS: "#BB1E10",
};
function hexToRgba(hex: string, alpha = 0.18): string {
  const h = hex.replace("#", "");
  const r = parseInt(h.slice(0, 2), 16);
  const g = parseInt(h.slice(2, 4), 16);
  const b = parseInt(h.slice(4, 6), 16);
  return `rgba(${r},${g},${b},${alpha})`;
}

type SeriesResp = {
  job_id: string;
  roots: string[];
  root?: string;
  series: { Parsimony: number[]; Dirichlet: number[]; HSS: number[] };
  counts: { Parsimony: number; Dirichlet: number; HSS: number };
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

  // Panel 1 (existing): pick one root; compare methods
  const [roots, setRoots] = useState<string[]>([]);
  const [root, setRoot] = useState<string>("");
  const [data, setData] = useState<SeriesResp | null>(null);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // Panel 2 (by-root for chosen method)
  const [methodSel, setMethodSel] = useState<MethodName>("Dirichlet");
  const [byRoot, setByRoot] = useState<Record<string, number[]> | null>(null);
  const [byRootLoading, setByRootLoading] = useState(false);
  const [byRootErr, setByRootErr] = useState<string | null>(null);

  // Simple in-memory cache: key = `${job}|final|${method}|${root}`
  const byRootCacheRef = useRef<Record<string, number[]>>({});

  // ---- EXPORT / DOWNLOAD STATE & HELPERS ----
  const exportPrefix = useMemo(() => buildExportPrefix?.({ job }) ?? (job || "violin"), [job]);

  const violinRefs = useRef<Record<string, PlotlyHTMLElement | null>>({});
  const registerViolin = useCallback(
    (id: string) =>
      (_figure: unknown, gd: PlotlyHTMLElement) => {
        violinRefs.current[id] = gd;
      },
    []
  );

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

  const utcStamp = () =>
    new Date().toISOString().replace(/\.\d{3}Z$/, "Z").replace(/[-:]/g, "").replace("T", "_");

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
    } catch { /* ignore */ }
  }, []);

  // ===== Panel 1: Fetch data (and roots) for current job/root =====
  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Please pick a job on the Precomputed Results page.");
      setData(null);
      setRoots([]);
      setRoot("");
      return;
    }
    setLoading(true);
    setErr(null);
    try {
      const u = new URL("/api/violin", window.location.origin);
      u.searchParams.set("job", job);
      u.searchParams.set("metric", "final");
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

      const resp = json as SeriesResp;

      setRoots(resp.roots ?? []);
      if ((!root || !resp.roots?.includes(root)) && (resp.root || resp.roots?.length)) {
        setRoot(resp.root ?? resp.roots[0]);
      }

      setData(resp);
    } catch (e: unknown) {
      setErr(e instanceof Error ? e.message : "Failed to load");
      setData(null);
      setRoots([]);
    } finally {
      setLoading(false);
    }
  }, [job, root]);

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

  // Build violin traces for panel 1 (root fixed; compare methods) — colored
  const violinTraces: Partial<Data>[] = useMemo(() => {
    if (!data) return [];
    return METHOD_ORDER.map((m) => {
      const color = COLORS[m];
      return {
        type: "violin",
        name: m,
        y: data.series[m],
        box: { visible: true },
        meanline: { visible: true },
        points: "all",
        jitter: 0.3,
        scalemode: "width",
        showlegend: false,
        line: { color },
        fillcolor: hexToRgba(color, 0.18),
        opacity: 0.95,
        marker: { color },
      } as Partial<Data>;
    });
  }, [data]);

  const medianLineTrace: Partial<Data> | null = useMemo(() => {
    if (!data) return null;
    const medians = METHOD_ORDER.map((m) => median(data.series[m]));
    return {
      type: "scatter",
      mode: "lines+markers",
      x: METHOD_ORDER as unknown as string[],
      y: medians,
      line: { color: "#111827", width: 2 }, // near-black, stands out from fills
      marker: { color: "#111827", size: 8 },
      name: "Median",
      showlegend: false,
      hoverinfo: "skip" as const,
    };
  }, [data]);

  const layout: Partial<Layout> = {
    title: { text: `ll_final · root=${root || "(none)"}`, font: { color: "black" } },
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

  // ====================== Panel 2 logic (by-root; metric=final) ======================

  const fetchSeriesForRoot = useCallback(
    async (rootName: string, method: MethodName): Promise<number[]> => {
      const key = `${job}|final|${method}|${rootName}`;
      const cached = byRootCacheRef.current[key];
      if (cached) return cached;

      const u = new URL("/api/violin", window.location.origin);
      u.searchParams.set("job", job);
      u.searchParams.set("metric", "final");
      u.searchParams.set("root", rootName);
      const res = await fetch(u.toString(), { cache: "no-store" });

      const ct = res.headers.get("content-type") || "";
      if (!ct.includes("application/json")) return [];
      const json: unknown = await res.json().catch(() => null);
      if (!res.ok || !json || isErr(json)) return [];
      const series = (json as SeriesResp).series as Record<MethodName, number[]>;
      const arr = Array.isArray(series?.[method]) ? series[method] : [];
      byRootCacheRef.current[key] = arr;
      return arr;
    },
    [job]
  );

  const loadByRoot = useCallback(async () => {
    if (!job || !roots.length) {
      setByRoot(null);
      return;
    }
    setByRootLoading(true);
    setByRootErr(null);
    try {
      const entries = await Promise.all(
        roots.map(async (r) => [r, await fetchSeriesForRoot(r, methodSel)] as const)
      );
      const obj: Record<string, number[]> = {};
      for (const [r, arr] of entries) {
        if (Array.isArray(arr) && arr.length) obj[r] = arr;
      }
      setByRoot(Object.keys(obj).length ? obj : null);
    } catch (e) {
      setByRootErr(e instanceof Error ? e.message : "Failed to load per-root series");
      setByRoot(null);
    } finally {
      setByRootLoading(false);
    }
  }, [job, roots, methodSel, fetchSeriesForRoot]);

  useEffect(() => {
    if (job) void loadByRoot();
  }, [job, roots, methodSel, loadByRoot]);

  // Colored by-root violins use the selected method's color
  const byRootTraces: Partial<Data>[] = useMemo(() => {
    if (!byRoot) return [];
    const color = COLORS[methodSel];
    const fill = hexToRgba(color, 0.18);
    return Object.entries(byRoot).map(([r, arr]) => ({
      type: "violin",
      name: r,
      y: arr,
      box: { visible: true },
      meanline: { visible: true },
      points: "all",
      jitter: 0.3,
      scalemode: "width",
      showlegend: false,
      line: { color },
      fillcolor: fill,
      opacity: 0.95,
      marker: { color },
    }));
  }, [byRoot, methodSel]);

  const layoutByRoot: Partial<Layout> = useMemo(
    () => ({
      title: { text: `ll_final · by root · method=${methodSel}`, font: { color: "black" } },
      xaxis: {
        title: { text: "Root", font: { color: "black" } },
        tickfont: { color: "black" },
        categoryorder: "array",
        categoryarray: roots,
        zerolinecolor: "#e5e7eb",
        gridcolor: "#f3f4f6",
      },
      yaxis: {
        title: { text: "ll_final", font: { color: "black" } },
        tickfont: { color: "black" },
        zerolinecolor: "#e5e7eb",
        gridcolor: "#f3f4f6",
      },
      margin: { l: 60, r: 20, t: 40, b: 80 },
      height: 520,
      showlegend: false,
      paper_bgcolor: "white",
      plot_bgcolor: "white",
    }),
    [methodSel, roots]
  );

  // ====================== RENDER ======================

  return (
    <div className="p-6 max-w-[1200px] mx-auto bg-white text-black">
      <h1 className="text-2xl font-bold mb-4">Violin Plots of Maximum Log Likelihood (ll_final)</h1>

      <div className="flex flex-wrap gap-3 items-end mb-4">
        <label className="flex flex-col">
          <span className="text-sm">Root node</span>
          <select
            className="border rounded px-3 py-2"
            value={root}
            onChange={(e) => setRoot(e.target.value)}
            disabled={!roots.length}
          >
            {roots.length ? (
              roots.map((n) => (
                <option key={n} value={n}>
                  {n}
                </option>
              ))
            ) : (
              <option value="" disabled>
                {job ? "No roots for this job" : "Select a job first"}
              </option>
            )}
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

      {/* Panel 1: one root, compare methods */}
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
                  <span
                    className="inline-block w-3 h-3 rounded-full"
                    style={{ backgroundColor: COLORS[m] }}
                    aria-hidden
                  />
                  <span>{m}</span>
                </div>
                <div className="text-sm">{data.counts[m]} repetitions</div>
              </div>
            ))}
          </div>

          <div className="bg-white border rounded p-2 mb-8">
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
          {job ? "No data for selected root." : "Select a job on the Precomputed Results page."}
        </div>
      )}

      {/* Panel 2 — by-root distribution for selected method */}
      <div className="flex flex-wrap gap-3 items-center mb-3">
        <label className="flex flex-col">
          <span className="text-sm">Method (by-root panel)</span>
          <select
            className="border rounded px-3 py-2"
            value={methodSel}
            onChange={(e) => setMethodSel(e.target.value as MethodName)}
            disabled={!roots.length}
          >
            {METHOD_ORDER.map((m) => (
              <option key={m} value={m}>{m}</option>
            ))}
          </select>
        </label>

        <div className="ml-auto flex gap-2">
          <button
            className="px-4 py-2 rounded bg-black text-white disabled:opacity-50"
            onClick={() => downloadViolin("ll_final_by_root", "png")}
            disabled={!byRoot || !Object.keys(byRoot).length}
            title="Download by-root violin as PNG"
          >
            Download PNG
          </button>
        </div>
      </div>

      {byRootLoading ? (
        <div className="p-4 border rounded">Loading by-root data…</div>
      ) : byRootErr ? (
        <div className="p-4 border rounded text-red-600">{byRootErr}</div>
      ) : byRoot && Object.keys(byRoot).length ? (
        <div className="bg-white border rounded p-2">
          <Plot
            data={byRootTraces}
            layout={layoutByRoot}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
            onInitialized={registerViolin("ll_final_by_root")}
            onUpdate={registerViolin("ll_final_by_root")}
          />
        </div>
      ) : (
        <div className="p-4 border rounded">No by-root data for the selected method.</div>
      )}
    </div>
  );
}
