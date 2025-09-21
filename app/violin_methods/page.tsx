// app/violin_methods/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import dynamic from "next/dynamic";
import type { Data, Layout, PlotlyHTMLElement } from "plotly.js";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

type MethodName = "Parsimony" | "Dirichlet" | "HSS";
const METHOD_ORDER: MethodName[] = ["Dirichlet", "Parsimony", "HSS"];

/** --- Shared color scheme --- */
const COLORS: Record<MethodName, string> = {
  Dirichlet: "#D2691E", // chocolate
  Parsimony: "#FF6B3D",
  HSS: "#BB1E10",
};
function hexToRgba(hex: string, alpha = 0.22) {
  const m = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  if (!m) return `rgba(0,0,0,${alpha})`;
  const r = parseInt(m[1], 16);
  const g = parseInt(m[2], 16);
  const b = parseInt(m[3], 16);
  return `rgba(${r}, ${g}, ${b}, ${alpha})`;
}
const showRoot = (s: string) => s.replace(/^h_/, "");

type SeriesResp = {
  job_id: string;
  roots: string[];
  root?: string;
  series: { Parsimony: number[]; Dirichlet: number[]; HSS: number[] };
  counts: { Parsimony: number; Dirichlet: number; HSS: number };
  metric: "final" | "init";
};
type ApiErr = { error: string };

function isErr(x: unknown): x is ApiErr {
  return !!x && typeof x === "object" && "error" in x;
}
function median(arr: number[]): number {
  const a = arr.filter((v) => Number.isFinite(v)).slice().sort((x, y) => x - y);
  const n = a.length;
  if (!n) return NaN;
  return n % 2 ? a[(n - 1) >> 1] : (a[n / 2 - 1] + a[n / 2]) / 2;
}

export default function ViolinMethodsPage() {
  const [job, setJob] = useState<string>("");
  const [roots, setRoots] = useState<string[]>([]);
  const [root, setRoot] = useState<string>("");

  // NEW: layer state + labels
  const [layer, setLayer] = useState<number>(1); // 0=coarse, 1=medium, 2=fine
  const LAYER_LABEL: Record<number, string> = { 0: "coarse (0)", 1: "medium (1)", 2: "fine (2)" };

  const [finalData, setFinalData] = useState<SeriesResp | null>(null);
  const [initData, setInitData] = useState<SeriesResp | null>(null);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // --- download helpers ---
  const plotsRef = useRef<Record<string, PlotlyHTMLElement | null>>({});
  const register = useCallback(
    (id: string) =>
      (_fig: unknown, gd: PlotlyHTMLElement) => {
        plotsRef.current[id] = gd;
      },
    []
  );
  const toImage = useCallback(async (gd: PlotlyHTMLElement) => {
    const w = window as unknown as { Plotly?: { toImage: (el: PlotlyHTMLElement, o: { format: "png"; scale: number }) => Promise<string> } };
    if (!w.Plotly?.toImage) throw new Error("Plotly.toImage not available");
    return w.Plotly.toImage(gd, { format: "png", scale: 2 });
  }, []);
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
  const downloadPNG = useCallback(
    async (id: "final_methods" | "init_methods") => {
      const gd = plotsRef.current[id];
      if (!gd) return;
      const url = await toImage(gd);
      const fname = `${job || "job"}__violin_methods__${id}__root-${root || "NA"}__layer-${layer}__t${utcStamp()}.png`;
      downloadURI(url, fname);
    },
    [job, root, layer, toImage]
  );

  // read job
  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId") ?? "";
      setJob(saved);
    } catch {}
  }, []);

  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Please pick a job on the Precomputed Results page.");
      setRoots([]); setRoot(""); setFinalData(null); setInitData(null);
      return;
    }
    setLoading(true); setErr(null);
    try {
      // one call to get roots & default root (final)
      const baseU = new URL("/api/violin", window.location.origin);
      baseU.searchParams.set("job", job);
      baseU.searchParams.set("metric", "final");
      baseU.searchParams.set("layer", String(layer)); // NEW
      const baseRes = await fetch(baseU.toString(), { cache: "no-store" });
      const baseJson: unknown = await baseRes.json();
      if (!baseRes.ok || isErr(baseJson)) throw new Error(isErr(baseJson) ? baseJson.error : `HTTP ${baseRes.status}`);
      const base = baseJson as SeriesResp;

      setRoots(base.roots ?? []);
      const chosenRoot = (root && base.roots.includes(root)) ? root : (base.root || base.roots[0] || "");
      setRoot(chosenRoot);

      // fetch both metrics for the selected root
      const [finalResp, initResp] = await Promise.all(
        (["final", "init"] as const).map(async (metric) => {
          const u = new URL("/api/violin", window.location.origin);
          u.searchParams.set("job", job);
          if (chosenRoot) u.searchParams.set("root", chosenRoot);
          u.searchParams.set("metric", metric);
          u.searchParams.set("layer", String(layer)); // NEW
          const r = await fetch(u.toString(), { cache: "no-store" });
          const j: unknown = await r.json();
          if (!r.ok || isErr(j)) throw new Error(isErr(j) ? j.error : `HTTP ${r.status}`);
          return j as SeriesResp;
        })
      );

      setFinalData(finalResp);
      setInitData(initResp);
    } catch (e: unknown) {
      setErr(e instanceof Error ? e.message : "Failed to load");
      setFinalData(null); setInitData(null); setRoots([]);
    } finally {
      setLoading(false);
    }
  }, [job, root, layer]); // NEW: layer dep

  useEffect(() => {
    if (job) void load();
  }, [job, root, layer, load]); // NEW: layer dep

  // traces builder with colors
  const buildMethodTraces = useCallback((d: SeriesResp | null): Partial<Data>[] => {
    if (!d) return [];
    return METHOD_ORDER.map((m) => ({
      type: "violin",
      name: m,
      y: d.series[m],
      box: { visible: true },
      meanline: { visible: true },
      points: "all",
      jitter: 0.3,
      scalemode: "width",
      showlegend: false,
      line: { color: COLORS[m] },
      fillcolor: hexToRgba(COLORS[m], 0.22),
      opacity: 0.95,
      marker: { color: COLORS[m] },
    }));
  }, []);

  const medianTrace = useCallback((d: SeriesResp | null): Partial<Data> | null => {
    if (!d) return null;
    const meds = METHOD_ORDER.map((m) => median(d.series[m]));
    return {
      type: "scatter",
      mode: "lines+markers",
      x: METHOD_ORDER as unknown as string[],
      y: meds,
      line: { color: "#111827", width: 2 },
      marker: { color: "#111827", size: 8 },
      name: "Median",
      showlegend: false,
      hoverinfo: "skip" as const,
    };
  }, []);

  const layout = useCallback((title: string): Partial<Layout> => ({
    title: { text: title, font: { color: "black" } },
    xaxis: {
      title: { text: "Method", font: { color: "black" } },
      tickfont: { color: "black" },
      categoryorder: "array",
      categoryarray: METHOD_ORDER,
      zerolinecolor: "#e5e7eb",
      gridcolor: "#f3f4f6",
    },
    yaxis: {
      title: { text: "log-likelihood", font: { color: "black" } },
      tickfont: { color: "black" },
      zerolinecolor: "#e5e7eb",
      gridcolor: "#f3f4f6",
    },
    margin: { l: 60, r: 20, t: 40, b: 50 },
    height: 420,
    showlegend: false,
    paper_bgcolor: "white",
    plot_bgcolor: "white",
  }), []);

  const titleFinal = useMemo(
    () => `Across Methods — final — root=${root ? showRoot(root) : "(none)"} — layer=${LAYER_LABEL[layer]}`,
    [root, layer]
  );
  const titleInit = useMemo(
    () => `Across Methods — init — root=${root ? showRoot(root) : "(none)"} — layer=${LAYER_LABEL[layer]}`,
    [root, layer]
  );

  return (
    <div className="p-6 max-w-[1200px] mx-auto bg-white text-black">
      <h1 className="text-2xl font-bold mb-4">Violin · Across Methods (init &amp; final)</h1>

      <div className="flex flex-wrap gap-3 items-end mb-4">
        <label className="flex flex-col">
          <span className="text-sm">Root</span>
          <select
            className="border rounded px-3 py-2"
            value={root}
            onChange={(e) => setRoot(e.target.value)}
            disabled={!roots.length}
          >
            {roots.length ? (
              roots.map((n) => <option key={n} value={n}>{n}</option>)
            ) : (
              <option value="" disabled>{job ? "No roots for this job" : "Select a job first"}</option>
            )}
          </select>
        </label>

        {/* NEW: Layer segmented control */}
        <div className="flex flex-col">
          <span className="text-sm">Layer</span>
          <div className="inline-flex rounded border overflow-hidden">
            {[0, 1, 2].map((lv) => (
              <button
                key={lv}
                type="button"
                onClick={() => setLayer(lv)}
                className={`px-3 py-2 text-sm border-r last:border-r-0 ${
                  layer === lv ? "bg-black text-white" : "bg-white text-black hover:bg-gray-100"
                }`}
                title={LAYER_LABEL[lv]}
              >
                {lv}
              </button>
            ))}
          </div>
        </div>

        <div className="text-sm ml-auto">
          {job ? <>Job: <span className="font-mono">{job}</span></> : "No job selected"}
          {loading && <span className="ml-3 animate-pulse">Loading…</span>}
          {err && <span className="ml-3 text-red-600">{err}</span>}
        </div>
      </div>

      {/* Final */}
      <section className="bg-white border rounded p-3 mb-6">
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg font-semibold">ll_final · root={root || "(none)"} · layer={layer}</h2>
          <button
            className="px-3 py-1 text-sm rounded bg-black text-white disabled:opacity-50"
            onClick={() => downloadPNG("final_methods")}
            disabled={!finalData}
            title="Download PNG"
          >
            Download PNG
          </button>
        </div>
        {finalData ? (
          <Plot
            data={[...buildMethodTraces(finalData), ...(medianTrace(finalData) ? [medianTrace(finalData)!] : [])]}
            layout={layout(titleFinal)}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
            onInitialized={register("final_methods")}
            onUpdate={register("final_methods")}
          />
        ) : (
          <div className="p-3 border rounded">No data</div>
        )}
      </section>

      {/* Init */}
      <section className="bg-white border rounded p-3">
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg font-semibold">ll_init · root={root || "(none)"} · layer={layer}</h2>
          <button
            className="px-3 py-1 text-sm rounded bg-black text-white disabled:opacity-50"
            onClick={() => downloadPNG("init_methods")}
            disabled={!initData}
            title="Download PNG"
          >
            Download PNG
          </button>
        </div>
        {initData ? (
          <Plot
            data={[...buildMethodTraces(initData), ...(medianTrace(initData) ? [medianTrace(initData)!] : [])]}
            layout={layout(titleInit)}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
            onInitialized={register("init_methods")}
            onUpdate={register("init_methods")}
          />
        ) : (
          <div className="p-3 border rounded">No data</div>
        )}
      </section>
    </div>
  );
}
