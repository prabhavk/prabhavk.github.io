// app/violin_roots/page.tsx
"use client";

import React, { useCallback, useEffect, useRef, useState } from "react";
import dynamic from "next/dynamic";
import type { Data, Layout, PlotlyHTMLElement } from "plotly.js";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

type MethodName = "Parsimony" | "Dirichlet" | "HSS";

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

export default function ViolinRootsPage() {
  const [job, setJob] = useState<string>("");
  const [roots, setRoots] = useState<string[]>([]);
  const [method, setMethod] = useState<MethodName>("Dirichlet");

  const [finalByRoot, setFinalByRoot] = useState<Record<string, number[]> | null>(null);
  const [initByRoot, setInitByRoot] = useState<Record<string, number[]> | null>(null);

  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // --- download helpers (one combined figure) ---
  const plotRef = useRef<PlotlyHTMLElement | null>(null);
  const register = useCallback((_fig: unknown, gd: PlotlyHTMLElement) => {
    plotRef.current = gd;
  }, []);
  const toImage = useCallback(async (gd: PlotlyHTMLElement) => {
    const w = window as unknown as {
      Plotly?: { toImage: (el: PlotlyHTMLElement, o: { format: "png"; scale: number }) => Promise<string> };
    };
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
  const downloadPNG = useCallback(async () => {
    const gd = plotRef.current;
    if (!gd) return;
    const url = await toImage(gd);
    const fname = `${job || "job"}__violin_roots__combined__method-${method}__t${utcStamp()}.png`;
    downloadURI(url, fname);
  }, [job, method, toImage]);

  // read job
  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId") ?? "";
      setJob(saved);
    } catch {}
  }, []);

  // load root list first (final metric is fine just to list roots)
  useEffect(() => {
    if (!job) {
      setRoots([]);
      return;
    }
    (async () => {
      try {
        const u = new URL("/api/violin", window.location.origin);
        u.searchParams.set("job", job);
        u.searchParams.set("metric", "final");
        const res = await fetch(u.toString(), { cache: "no-store" });
        const j: unknown = await res.json();
        if (!res.ok || isErr(j)) throw new Error(isErr(j) ? j.error : `HTTP ${res.status}`);
        const resp = j as SeriesResp;
        setRoots(resp.roots ?? []);
      } catch (e: any) {
        setErr(e?.message ?? "Failed to load roots");
        setRoots([]);
      }
    })();
  }, [job]);

  // fetch arrays per root for a metric & method
  const fetchMetricByRoot = useCallback(
    async (metric: "final" | "init", m: MethodName) => {
      const out: Record<string, number[]> = {};
      for (const r of roots) {
        const u = new URL("/api/violin", window.location.origin);
        u.searchParams.set("job", job);
        u.searchParams.set("root", r);
        u.searchParams.set("metric", metric);
        const res = await fetch(u.toString(), { cache: "no-store" });
        const j: unknown = await res.json();
        if (!res.ok || isErr(j)) continue;
        const series = (j as SeriesResp).series as Record<MethodName, number[]>;
        out[r] = Array.isArray(series[m]) ? series[m] : [];
      }
      const filtered: Record<string, number[]> = {};
      for (const [r, arr] of Object.entries(out)) {
        if (arr && arr.length) filtered[r] = arr;
      }
      return filtered;
    },
    [job, roots]
  );

  // load both metrics when method or roots change
  useEffect(() => {
    if (!job || !roots.length) {
      setFinalByRoot(null);
      setInitByRoot(null);
      return;
    }
    (async () => {
      setLoading(true);
      setErr(null);
      try {
        const [f, i] = await Promise.all([fetchMetricByRoot("final", method), fetchMetricByRoot("init", method)]);
        setFinalByRoot(Object.keys(f).length ? f : null);
        setInitByRoot(Object.keys(i).length ? i : null);
      } catch (e: any) {
        setErr(e?.message ?? "Failed to load by-root series");
        setFinalByRoot(null);
        setInitByRoot(null);
      } finally {
        setLoading(false);
      }
    })();
  }, [job, roots, method, fetchMetricByRoot]);

  /** Build combined side-by-side violins (init left, final right) + median connectors
   *  Adds replicate info via `customdata` and `hovertemplate`.
   */
  const buildCombinedTraces = useCallback((): Partial<Data>[] => {
    if (!finalByRoot && !initByRoot) return [];

    // Use the order from `roots`, but only include those with at least one series
    const ordered = roots.filter((r) => finalByRoot?.[r]?.length || initByRoot?.[r]?.length);

    const c = COLORS[method];
    const fillInit = hexToRgba(c, 0.16);
    const fillFinal = hexToRgba(c, 0.28);

    // numeric x positions with a small offset for side-by-side placement
    const d = 0.18; // horizontal offset

    const traces: Partial<Data>[] = [];

    // 1) Violin traces (init left, final right) — include replicate index
    ordered.forEach((r, idx) => {
      const i = idx + 1; // center position for this root
      const arrInit = initByRoot?.[r] ?? [];
      const arrFinal = finalByRoot?.[r] ?? [];

      if (arrInit.length) {
        const reps = arrInit.map((_, k) => k + 1); // 1-based replicate index
        traces.push({
          type: "violin",
          name: `${r} (init)`,
          y: arrInit,
          x: new Array(arrInit.length).fill(i - d),
          customdata: reps,
          hovertemplate: `root=${r}<br>metric=init<br>rep=%{customdata}<br>ll=%{y}<extra></extra>`,
          box: { visible: true }, // median line inside the box
          meanline: { visible: true },
          points: "all",
          jitter: 0.3,
          scalemode: "width",
          showlegend: false,
          line: { color: c },
          fillcolor: fillInit,
          opacity: 0.9,
          marker: { color: c },
        } as Partial<Data>);
      }

      if (arrFinal.length) {
        const reps = arrFinal.map((_, k) => k + 1); // 1-based replicate index
        traces.push({
          type: "violin",
          name: `${r} (final)`,
          y: arrFinal,
          x: new Array(arrFinal.length).fill(i + d),
          customdata: reps,
          hovertemplate: `root=${r}<br>metric=final<br>rep=%{customdata}<br>ll=%{y}<extra></extra>`,
          box: { visible: true }, // median line inside the box
          meanline: { visible: true },
          points: "all",
          jitter: 0.3,
          scalemode: "width",
          showlegend: false,
          line: { color: c },
          fillcolor: fillFinal,
          opacity: 0.95,
          marker: { color: c },
        } as Partial<Data>);
      }
    });

    // 2) Median connectors (short segment from init median to final median)
    const connX: (number | null)[] = [];
    const connY: (number | null)[] = [];

    ordered.forEach((r, idx) => {
      const i = idx + 1;
      const arrInit = initByRoot?.[r] ?? [];
      const arrFinal = finalByRoot?.[r] ?? [];
      if (!arrInit.length || !arrFinal.length) return;

      const mInit = median(arrInit);
      const mFinal = median(arrFinal);

      // segment (i - d, mInit) -> (i + d, mFinal), then a null to break
      connX.push(i - d, i + d, null);
      connY.push(mInit, mFinal, null);
    });

    if (connX.length) {
      traces.push({
        type: "scatter",
        mode: "lines+markers",
        x: connX,
        y: connY,
        line: { color: "black", width: 2 },
        marker: { color: "black", size: 6 },
        hoverinfo: "skip",
        showlegend: false,
      } as Partial<Data>);
    }

    return traces;
  }, [finalByRoot, initByRoot, method, roots]);

  const layout: Partial<Layout> = {
    title: { text: `Across Roots — init (left) vs final (right) · ${method}`, font: { color: "black" } },
    xaxis: {
      title: { text: "Root", font: { color: "black" } },
      tickfont: { color: "black" },
      type: "linear",
      tickmode: "array",
      tickvals: roots.map((_, idx) => idx + 1),
      ticktext: roots,
      zerolinecolor: "#e5e7eb",
      gridcolor: "#f3f4f6",
    },
    yaxis: {
      title: { text: "log-likelihood", font: { color: "black" } },
      tickfont: { color: "black" },
      zerolinecolor: "#e5e7eb",
      gridcolor: "#f3f4f6",
    },
    margin: { l: 60, r: 20, t: 40, b: 80 },
    height: 560,
    showlegend: false,
    paper_bgcolor: "white",
    plot_bgcolor: "white",
    hovermode: "closest",
  };

  return (
    <div className="p-6 max-w-[1200px] mx-auto bg-white text-black">
      <h1 className="text-2xl font-bold mb-4">Violin · Across Roots (init &amp; final on one plot)</h1>

      <div className="flex flex-wrap gap-3 items-end mb-4">
        <label className="flex flex-col">
          <span className="text-sm">Method</span>
          <select
            className="border rounded px-3 py-2"
            value={method}
            onChange={(e) => setMethod(e.target.value as MethodName)}
            disabled={!roots.length}
          >
            {(["Dirichlet", "Parsimony", "HSS"] as MethodName[]).map((m) => (
              <option key={m} value={m}>
                {m}
              </option>
            ))}
          </select>
        </label>

        <div className="text-sm ml-auto">
          {job ? (
            <>
              Job: <span className="font-mono">{job}</span>
            </>
          ) : (
            "No job selected"
          )}
          {loading && <span className="ml-3 animate-pulse">Loading…</span>}
          {err && <span className="ml-3 text-red-600">{err}</span>}
        </div>
      </div>

      <section className="bg-white border rounded p-3">
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg font-semibold">Init (left) vs Final (right)</h2>
          <button
            className="px-3 py-1 text-sm rounded bg-black text-white disabled:opacity-50"
            onClick={downloadPNG}
            disabled={!finalByRoot && !initByRoot}
            title="Download PNG"
          >
            Download PNG
          </button>
        </div>

        {!finalByRoot && !initByRoot ? (
          <div className="p-3 border rounded">No data</div>
        ) : (
          <Plot
            data={buildCombinedTraces()}
            layout={layout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
            onInitialized={register}
            onUpdate={register}
          />
        )}
      </section>
    </div>
  );
}
