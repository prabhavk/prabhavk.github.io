// app/ecdll/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import dynamic from "next/dynamic";
import type { Layout, PlotData } from "plotly.js";
import type { PlotParams } from "react-plotly.js";

const Plot = dynamic<PlotParams>(() => import("react-plotly.js"), { ssr: false });

/* ---------------- Types ---------------- */

type MethodKey = "dirichlet" | "parsimony" | "hss";
type LayerIdx = 0 | 1 | 2;

const METHOD_LABEL: Record<MethodKey, string> = {
  dirichlet: "Dirichlet",
  parsimony: "Parsimony",
  hss: "HSS",
};

const LAYER_LABEL: Record<LayerIdx, string> = {
  0: "Bottom",
  1: "Medium",
  2: "Top",
};

const COLORS: Record<MethodKey, string> = {
  dirichlet: "#D2691E", // chocolate
  parsimony: "#FF6B3D",
  hss: "#BB1E10",
};

type SummaryRow = {
  root: string | null;
  num_iterations: number;
  ll_init: number | null;
  ecd_ll_first: number | null;
  ecd_ll_final: number | null;
  ll_final: number | null;
  root_prob_final: number[] | null;
};

type Pt = [number, number]; // [iter, ecd_ll]

type ApiOk = {
  ok: true;
  job_id: string;
  method: MethodKey;
  layer: LayerIdx;
  rep: number | null;
  points: Pt[];
  summary: SummaryRow | null;
};

type ApiErr = { ok: false; error: string };

type ViolinRootsResp =
  | {
      job_id: string;
      roots: string[];
      root?: string;
      metric: "final" | "init";
    }
  | { error: string };

/* ---------------- Helpers ---------------- */

function normalizePoints(x: unknown): Pt[] {
  if (!Array.isArray(x)) return [];
  const out: Pt[] = [];
  for (const row of x) {
    if (Array.isArray(row) && row.length >= 2) {
      const it = Number(row[0]);
      const ll = Number(row[1]);
      if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
    } else if (row && typeof row === "object") {
      const it = Number((row as Record<string, unknown>).iter);
      const ll = Number((row as Record<string, unknown>).ll);
      if (Number.isFinite(it) && Number.isFinite(ll)) out.push([Math.round(it), ll]);
    }
  }
  out.sort((a, b) => a[0] - b[0]);
  return out;
}

/** Gap^α transform for ECD-LL curves (composite plot only). */
function transformGapPower(points: Pt[], alpha: number = 0.5, eps = 1e-9): Pt[] {
  if (!Array.isArray(points) || points.length === 0) return points;
  const ys = points.map(([, ll]) => ll).filter((v) => Number.isFinite(v));
  if (ys.length === 0) return points;

  const best = Math.max(...ys);
  const gaps = points.map(([, ll]) => best - ll);
  const gMax = Math.max(...gaps.filter((g) => Number.isFinite(g)), 0);

  return points.map(([it], i) => {
    const g = gaps[i];
    if (!Number.isFinite(g)) return [it, 0];
    const y = Math.pow(g / (gMax + eps), alpha);
    return [it, y];
  });
}

/** Solid line + markers (so single-point series are visible) */
function monoTrace(name: string, pts: Pt[] | null, color: string): Partial<PlotData> | null {
  if (!pts || !pts.length) return null;
  return {
    type: "scatter",
    mode: "lines+markers",
    name,
    x: pts.map((p) => p[0]),
    y: pts.map((p) => p[1]),
    line: { color, width: 2 },
    marker: { color, size: 5, line: { color, width: 0 } },
  };
}

/* ---------------- Page ---------------- */

export default function EcdllPage() {
  const [job, setJob] = useState<string>("");

  // NEW: layer and selector
  const [layer, setLayer] = useState<LayerIdx>(0);

  // roots and selected root
  const [roots, setRoots] = useState<string[]>([]);
  const [root, setRoot] = useState<string>("");

  const [err, setErr] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  const [reps, setReps] = useState<number[]>([]);
  const [rep, setRep] = useState<number | null>(null);

  const [dirPts, setDirPts] = useState<Pt[] | null>(null);
  const [parPts, setParPts] = useState<Pt[] | null>(null);
  const [hssPts, setSshPts] = useState<Pt[] | null>(null);

  // Now *used* in the summary table
  const [dirSum, setDirSum] = useState<SummaryRow | null>(null);
  const [parSum, setParSum] = useState<SummaryRow | null>(null);
  const [hssSum, setSshSum] = useState<SummaryRow | null>(null);

  const [alpha, setAlpha] = useState<number>(0.5);

  const graphRefs = useRef<Record<string, HTMLDivElement | null>>({
    composite: null,
    all: null,   // overlay plot ref
    dir: null,
    par: null,
    hss: null,
  });

  // Card/plot frame classes
  const cardCls = "border-2 border-black rounded p-3";
  const plotFrameCls = "border-2 border-black rounded";

  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId");
      if (saved) setJob(saved);
      const savedLayer = localStorage.getItem("emtr:ecdll:selectedLayer");
      if (savedLayer != null) {
        const n = Number(savedLayer);
        if (n === 0 || n === 1 || n === 2) setLayer(n as LayerIdx);
      }
    } catch {
      /* ignore */
    }
  }, []);

  useEffect(() => {
    try {
      localStorage.setItem("emtr:ecdll:selectedLayer", String(layer));
    } catch { /* ignore */ }
  }, [layer]);

  // Load roots for this job (reuse /api/violin list endpoint)
  const loadRoots = useCallback(async (jobId: string) => {
    const u = new URL("/api/violin", window.location.origin);
    u.searchParams.set("job", jobId);
    u.searchParams.set("metric", "final");
    const res = await fetch(u.toString(), { cache: "no-store" });
    const ct = (res.headers.get("content-type") || "").toLowerCase();
    const body = await res.text();
    if (!ct.includes("application/json")) {
      throw new Error(`Expected JSON for roots, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`);
    }
    const j = JSON.parse(body) as ViolinRootsResp;
    if ("error" in j) throw new Error(j.error);
    const list = Array.isArray(j.roots) ? j.roots : [];
    setRoots(list);
    // prefer backend-chosen default; else first available
    const preferred = j.root && list.includes(j.root) ? j.root : list[0] ?? "";
    setRoot((prev) => (prev && list.includes(prev) ? prev : preferred));
  }, []);

  // Load reps; pass root and layer (backend may ignore layer if not implemented)
  const loadReps = useCallback(
    async (jobId: string, rootSel?: string, layerSel?: LayerIdx) => {
      const u = new URL("/api/ecdll/reps", window.location.origin);
      u.searchParams.set("job_id", jobId);
      if (rootSel) u.searchParams.set("root", rootSel);
      if (layerSel !== undefined) u.searchParams.set("layer", String(layerSel));
      const res = await fetch(u.toString(), { cache: "no-store" });
      const j = (await res.json()) as { ok: boolean; reps?: number[]; error?: string };
      if (!res.ok || !j.ok) throw new Error(j.error || `HTTP ${res.status}`);
      const list = j.reps ?? [];
      setReps(list);
      if (list.length && (rep === null || !list.includes(rep))) setRep(list[0]);
    },
    [rep]
  );

  // Kick off roots + reps when job OR layer changes
  useEffect(() => {
    if (!job) return;
    (async () => {
      setErr(null);
      try {
        await loadRoots(job);
      } catch (e) {
        setErr(e instanceof Error ? e.message : "Failed to load roots");
        setRoots([]);
        setRoot("");
      }
      try {
        if (root) await loadReps(job, root, layer);
      } catch (e) {
        setErr(e instanceof Error ? e.message : "Failed to load repetitions");
        setReps([]);
        setRep(null);
      }
    })();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [job, layer, loadRoots, loadReps]);

  // When root changes, (re)load reps for this layer
  useEffect(() => {
    if (!job || !root) return;
    (async () => {
      try {
        await loadReps(job, root, layer);
      } catch (e) {
        setErr(e instanceof Error ? e.message : "Failed to load repetitions");
        setReps([]);
        setRep(null);
      }
    })();
  }, [job, root, layer, loadReps]);

  const fetchOne = useCallback(
    async (method: MethodKey, repSel: number | null, rootSel: string, layerSel: LayerIdx): Promise<{ pts: Pt[]; sum: SummaryRow | null }> => {
      const u = new URL("/api/ecdll", window.location.origin);
      u.searchParams.set("job_id", job);
      u.searchParams.set("method", method);
      u.searchParams.set("layer", String(layerSel)); // NEW
      if (repSel !== null) u.searchParams.set("rep", String(repSel));
      if (rootSel) u.searchParams.set("root", rootSel);
      const res = await fetch(u.toString(), { cache: "no-store" });
      const ct = res.headers.get("content-type") || "";
      if (!ct.toLowerCase().includes("application/json")) {
        const body = await res.text().catch(() => "");
        throw new Error(
          `Expected JSON for ${method}, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`
        );
      }
      const j = (await res.json()) as ApiOk | ApiErr;

      if (!res.ok || !("ok" in j) || j.ok === false) {
        const msg = "error" in (j as ApiErr) ? (j as ApiErr).error : `HTTP ${res.status}`;
        throw new Error(`${METHOD_LABEL[method]}: ${msg}`);
      }
      const ok = j as ApiOk;
      return { pts: normalizePoints(ok.points), sum: ok.summary };
    },
    [job]
  );

  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Choose a job on the Precomputed Results page.");
      setDirPts(null);
      setParPts(null);
      setSshPts(null);
      setDirSum(null);
      setParSum(null);
      setSshSum(null);
      return;
    }
    if (rep === null) {
      setErr("No repetition selected.");
      return;
    }
    if (!root) {
      setErr("No root selected.");
      return;
    }
    setErr(null);
    setLoading(true);
    try {
      const [d, p, s] = await Promise.allSettled([
        fetchOne("dirichlet", rep, root, layer),
        fetchOne("parsimony", rep, root, layer),
        fetchOne("hss", rep, root, layer),
      ]);

      setDirPts(d.status === "fulfilled" ? d.value.pts : []);
      setParPts(p.status === "fulfilled" ? p.value.pts : []);
      setSshPts(s.status === "fulfilled" ? s.value.pts : []);

      setDirSum(d.status === "fulfilled" ? d.value.sum : null);
      setParSum(p.status === "fulfilled" ? p.value.sum : null);
      setSshSum(s.status === "fulfilled" ? s.value.sum : null);

      const errs: string[] = [];
      if (d.status === "rejected") errs.push((d as PromiseRejectedResult).reason?.message || "Dirichlet fetch failed");
      if (p.status === "rejected") errs.push((p as PromiseRejectedResult).reason?.message || "Parsimony fetch failed");
      if (s.status === "rejected") errs.push((s as PromiseRejectedResult).reason?.message || "HSS fetch failed");
      if (errs.length) setErr(errs.join(" | "));
    } catch (e) {
      setErr(e instanceof Error ? e.message : "Failed to load ECD-LL data");
      setDirPts([]);
      setParPts([]);
      setSshPts([]);
      setDirSum(null);
      setParSum(null);
      setSshSum(null);
    } finally {
      setLoading(false);
    }
  }, [job, rep, root, layer, fetchOne]);

  useEffect(() => {
    if (job && rep !== null && root) void load();
  }, [job, rep, root, layer, load]);

  /* -------- Layouts (integer x ticks everywhere) -------- */

  const baseLayout: Partial<Layout> = useMemo(
    () => ({
      margin: { l: 60, r: 20, t: 30, b: 50 },
      xaxis: {
        title: { text: "iteration" },
        tickmode: "linear",
        dtick: 1,
        tickformat: "d",
        ticks: "outside",
        showline: true,
        linecolor: "#000000",
        linewidth: 2,
        mirror: true,
        zeroline: false,
      },
      yaxis: {
        title: { text: "ECD-LL" },
        ticks: "outside",
        showline: true,
        linecolor: "#000000",
        linewidth: 2,
        mirror: true,
        zeroline: false,
      },
      height: 320,
      showlegend: true,
      plot_bgcolor: "white",
      paper_bgcolor: "white",
    }),
    []
  );

  const compositeLayout: Partial<Layout> = useMemo(
    () => ({
      ...baseLayout,
      yaxis: { title: { text: `Gap^α (0 = best)` } },
      showlegend: true,
      height: 360,
    }),
    [baseLayout]
  );

  /* -------- Build traces -------- */

  // Gap^α composite
  const compositeTraces = useMemo(() => {
    const traces: Partial<PlotData>[] = [];
    const push = (name: string, pts: Pt[] | null, color: string) => {
      if (!pts || !pts.length) return;
      const t = transformGapPower(pts, alpha);
      traces.push({
        type: "scatter",
        mode: "lines+markers",
        name,
        x: t.map((p) => p[0]),
        y: t.map((p) => p[1]),
        line: { color, width: 2 },
        marker: { color, size: 5, line: { color, width: 0 } },
      });
    };
    push("Dirichlet", dirPts, COLORS.dirichlet);
    push("Parsimony", parPts, COLORS.parsimony);
    push("HSS", hssPts, COLORS.hss);
    return traces;
  }, [dirPts, parPts, hssPts, alpha]);

  // Individual
  const dirTrace = useMemo(() => monoTrace("Dirichlet", dirPts, COLORS.dirichlet), [dirPts]);
  const parTrace = useMemo(() => monoTrace("Parsimony", parPts, COLORS.parsimony), [parPts]);
  const hssTrace = useMemo(() => monoTrace("HSS", hssPts, COLORS.hss), [hssPts]);

  // Overlay all three methods (raw ECD-LL)
  const allMethodsTraces = useMemo(() => {
    const arr: Partial<PlotData>[] = [];
    const a = monoTrace("Dirichlet", dirPts, COLORS.dirichlet);
    const b = monoTrace("Parsimony", parPts, COLORS.parsimony);
    const c = monoTrace("HSS", hssPts, COLORS.hss);
    if (a) arr.push(a);
    if (b) arr.push(b);
    if (c) arr.push(c);
    return arr;
  }, [dirPts, parPts, hssPts]);

  /* -------- Summary table data -------- */

  const tableRows = useMemo(() => {
    const rows: Array<{
      method: string;
      color: string;
      sum: SummaryRow | null;
    }> = [
      { method: "Dirichlet", color: COLORS.dirichlet, sum: dirSum },
      { method: "Parsimony", color: COLORS.parsimony, sum: parSum },
      { method: "HSS",       color: COLORS.hss,       sum: hssSum },
    ];
    return rows;
  }, [dirSum, parSum, hssSum]);

  /* -------- Download helpers -------- */

  type PlotKey = "composite" | "all" | "dir" | "par" | "hss";

  interface PlotlyLike {
    downloadImage: (
      gd: HTMLElement,
      opts: { format?: "png" | "jpeg" | "svg" | "webp"; width?: number; height?: number; filename?: string }
    ) => Promise<string>;
  }

  function isRecord(v: unknown): v is Record<string, unknown> {
    return typeof v === "object" && v !== null;
  }
  function isPlotlyLike(v: unknown): v is PlotlyLike {
    return isRecord(v) && typeof v["downloadImage"] === "function";
  }

  const resolvePlotlyFrom = useCallback(async (gd: HTMLDivElement): Promise<PlotlyLike> => {
    const win = (gd.ownerDocument?.defaultView ?? null) as (Window & { Plotly?: unknown }) | null;
    if (win?.Plotly && isPlotlyLike(win.Plotly)) return win.Plotly;
    const mod = (await import("plotly.js-dist")) as unknown;
    const candidate = isRecord(mod) && "default" in mod ? ((mod as { default: unknown }).default as unknown) : mod;
    if (isPlotlyLike(candidate)) return candidate;
    throw new Error("Plotly not found. Install 'plotly.js-dist' or expose window.Plotly.");
  }, []);

  const downloadFigure = useCallback(
    async (key: PlotKey, filename: string, w = 1600, h = 900) => {
      const gd = graphRefs.current[key];
      if (!gd) return;
      const Plotly = await resolvePlotlyFrom(gd);
      await Plotly.downloadImage(gd, { format: "png", width: w, height: h, filename });
    },
    [resolvePlotlyFrom]
  );

  const downloadAll = useCallback(async () => {
    const tasks: Array<Promise<void>> = [];
    const base = `ecdll_job-${job || "unknown"}_root-${(root || "na").replace(/[^A-Za-z0-9_-]+/g, "_")}_rep-${
      rep ?? "na"
    }_layer-${layer}`;

    if (graphRefs.current.composite)
      tasks.push(downloadFigure("composite", `${base}_composite_alpha-${alpha.toFixed(2)}`));
    if (graphRefs.current.all) tasks.push(downloadFigure("all", `${base}_all-methods`));
    if (graphRefs.current.dir) tasks.push(downloadFigure("dir", `${base}_dirichlet`));
    if (graphRefs.current.par) tasks.push(downloadFigure("par", `${base}_parsimony`));
    if (graphRefs.current.hss) tasks.push(downloadFigure("hss", `${base}_hss`));

    for (const t of tasks) await t;
  }, [alpha, job, root, rep, layer, downloadFigure]);

  /* ---------------- Render ---------------- */

  return (
    <div className="p-6 max-w-[1400px] mx-auto space-y-5">
      <div className="flex flex-wrap items-end gap-3">
        <h1 className="text-2xl text-black font-bold">Expected Complete-Data LL Convergence</h1>

        <div className="text-sm text-black ml-2">
          Job: <span className="font-mono">{job || "(none)"} </span>
        </div>

        {/* NEW: Layer selector */}
        <div className="flex items-center gap-2 ml-2">
          <span className="text-sm text-black">Layer</span>
          <div className="inline-flex rounded-md border overflow-hidden">
            {([0, 1, 2] as LayerIdx[]).map((L) => (
              <button
                key={L}
                type="button"
                onClick={() => setLayer(L)}
                className={`px-3 py-1.5 text-sm ${
                  layer === L ? "bg-black text-white" : "bg-white text-black"
                } ${L !== 2 ? "border-r" : ""}`}
                title={LAYER_LABEL[L]}
              >
                {LAYER_LABEL[L]}
              </button>
            ))}
          </div>
        </div>

        {/* Root selector */}
        <label className="ml-4 text-sm flex text-black items-center gap-2">
          Root:
          <select className="border rounded px-2 py-1 text-sm" value={root} onChange={(e) => setRoot(e.target.value)}>
            {roots.map((r) => (
              <option key={r} value={r}>
                {r.replace(/^h_/, "")}
              </option>
            ))}
          </select>
        </label>

        {/* Repetition selector */}
        <label className="text-sm flex text-black items-center gap-2">
          Rep:
          <select
            className="border rounded px-2 py-1 text-sm"
            value={rep ?? ""}
            onChange={(e) => setRep(e.target.value ? Number(e.target.value) : null)}
          >
            {reps.map((r) => (
              <option key={r} value={r}>
                {r}
              </option>
            ))}
          </select>
        </label>

        <button
          type="button"
          onClick={() => void load()}
          className="ml-auto px-3 py-1.5 rounded bg-black text-white disabled:opacity-50 text-sm"
          disabled={!job || rep === null || !root || loading || reps.length === 0}
        >
          {loading ? "Loading…" : "Reload"}
        </button>

        {/* Download ALL */}
        <button
          type="button"
          onClick={() => void downloadAll()}
          className="px-3 py-1.5 rounded bg-black text-white disabled:opacity-50 text-sm"
          disabled={loading}
          title="Download composite, overlay, and individuals as PNGs"
        >
          Download all PNGs
        </button>
      </div>

      {/* Alpha control for composite transform */}
      <div className="flex items-center gap-3">
        <label className="text-black text-sm">Composite transform α (0.2–0.9)</label>
        <input
          type="range"
          min={0.2}
          max={0.9}
          step={0.05}
          value={alpha}
          onChange={(e) => setAlpha(Number(e.target.value))}
        />
        <div className="text-xs text-gray-500">α = {alpha.toFixed(2)}</div>
      </div>

      {/* Composite plot (Gap^α) */}
      <section className={cardCls}>
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg text-black font-semibold">
            Composite (Gap^α transform)
            {root ? ` — root=${root.replace(/^h_/, "")}` : ""} · layer={LAYER_LABEL[layer]}
          </h2>
          <button
            type="button"
            onClick={() =>
              downloadFigure(
                "composite",
                `ecdll_job-${job || "unknown"}_root-${(root || "na").replace(/[^A-Za-z0-9_-]+/g, "_")}_rep-${
                  rep ?? "na"
                }_layer-${layer}_composite_alpha-${alpha.toFixed(2)}`
              )
            }
            className="px-2 py-1 rounded bg-black text-white"
          >
            Download PNG
          </button>
        </div>
        {dirPts?.length || parPts?.length || hssPts?.length ? (
          <div className={plotFrameCls}>
            <Plot
              data={compositeTraces}
              layout={compositeLayout}
              config={{ displayModeBar: false, responsive: true }}
              style={{ width: "100%", height: "100%" }}
              onInitialized={(_, gd) => {
                graphRefs.current.composite = gd as HTMLDivElement;
              }}
              onUpdate={(_, gd) => {
                graphRefs.current.composite = gd as HTMLDivElement;
              }}
            />
          </div>
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data found for any method.</div>
        )}
      </section>

      {/* Overlay all three methods (raw ECD-LL) */}
      <section className={cardCls}>
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg text-black font-semibold">
            All methods (ECD-LL)
            {root ? ` — root=${root.replace(/^h_/, "")}` : ""} · layer={LAYER_LABEL[layer]}
          </h2>
          <button
            type="button"
            onClick={() =>
              downloadFigure(
                "all",
                `ecdll_job-${job || "unknown"}_root-${(root || "na").replace(/[^A-Za-z0-9_-]+/g, "_")}_rep-${
                  rep ?? "na"
                }_layer-${layer}_all-methods`
              )
            }
            className="px-2 py-1 rounded bg-black text-white"
          >
            Download PNG
          </button>
        </div>
        {allMethodsTraces.length ? (
          <div className={plotFrameCls}>
            <Plot
              data={allMethodsTraces as PlotData[]}
              layout={baseLayout}
              config={{ displayModeBar: false, responsive: true }}
              style={{ width: "100%", height: "100%" }}
              onInitialized={(_, gd) => {
                graphRefs.current.all = gd as HTMLDivElement;
              }}
              onUpdate={(_, gd) => {
                graphRefs.current.all = gd as HTMLDivElement;
              }}
            />
          </div>
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data to overlay.</div>
        )}
      </section>

      {/* Summary table */}
      <section className={`${cardCls} text-black`}>
        <h2 className="text-lg font-semibold mb-2">
          Summary — rep={rep ?? "?"}{root ? ` · root=${root.replace(/^h_/, "")}` : ""} · layer={LAYER_LABEL[layer]}
        </h2>
        <div className="overflow-x-auto">
          <table className="min-w-full text-sm border-collapse text-black">
            <thead className="text-black">
              <tr className="border-b">
                <th className="text-left py-2 pr-3 font-medium">Method</th>
                <th className="text-left py-2 pr-3 font-medium"># iters</th>
                <th className="text-left py-2 pr-3 font-medium">init_ll</th>
                <th className="text-left py-2 pr-3 font-medium">ecd_ll_first</th>
                <th className="text-left py-2 pr-3 font-medium">ecd_ll_final</th>
                <th className="text-left py-2 pr-3 font-medium">final_ll</th>
                <th className="text-left py-2 pr-3 font-medium">root_prob_final [4]</th>
              </tr>
            </thead>
            <tbody className="text-black">
              {[
                { method: "Dirichlet", color: COLORS.dirichlet, sum: dirSum },
                { method: "Parsimony", color: COLORS.parsimony, sum: parSum },
                { method: "HSS",       color: COLORS.hss,       sum: hssSum },
              ].map(({ method, color, sum }) => (
                <tr key={method} className="border-b align-top">
                  <td className="py-2 pr-3">
                    <span className="inline-flex items-center gap-2">
                      <span className="inline-block w-3 h-3 rounded-sm" style={{ background: color }} />
                      {method}
                    </span>
                  </td>
                  {sum ? (
                    <>
                      <td className="py-2 pr-3">{sum.num_iterations ?? 0}</td>
                      <td className="py-2 pr-3">{sum.ll_init != null ? sum.ll_init.toFixed(6) : "—"}</td>
                      <td className="py-2 pr-3">{sum.ecd_ll_first != null ? sum.ecd_ll_first.toFixed(6) : "—"}</td>
                      <td className="py-2 pr-3">{sum.ecd_ll_final != null ? sum.ecd_ll_final.toFixed(6) : "—"}</td>
                      <td className="py-2 pr-3">{sum.ll_final != null ? sum.ll_final.toFixed(6) : "—"}</td>
                      <td className="py-2 pr-3 font-mono">
                        {Array.isArray(sum.root_prob_final) && sum.root_prob_final.length === 4
                          ? `[${sum.root_prob_final.map((p) => (Number.isFinite(p) ? p.toFixed(4) : "NaN")).join(", ")}]`
                          : "—"}
                      </td>
                    </>
                  ) : (
                    <td className="py-2 pr-3" colSpan={6}>No data</td>
                  )}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </section>

      {/* Individual plots */}
      <section className={cardCls}>
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg text-black font-semibold">
            Dirichlet{root ? ` — root=${root.replace(/^h_/, "")}` : ""} · layer={LAYER_LABEL[layer]}
          </h2>
          <button
            type="button"
            onClick={() =>
              downloadFigure(
                "dir",
                `ecdll_job-${job || "unknown"}_root-${(root || "na").replace(/[^A-Za-z0-9_-]+/g, "_")}_rep-${
                  rep ?? "na"
                }_layer-${layer}_dirichlet`
              )
            }
            className="px-2 py-1 rounded bg-black text-white"
          >
            Download PNG
          </button>
        </div>
        {dirTrace ? (
          <div className={plotFrameCls}>
            <Plot
              data={[dirTrace as Partial<PlotData>]}
              layout={baseLayout}
              config={{ displayModeBar: false, responsive: true }}
              style={{ width: "100%", height: "100%" }}
              onInitialized={(_, gd) => {
                graphRefs.current.dir = gd as HTMLDivElement;
              }}
              onUpdate={(_, gd) => {
                graphRefs.current.dir = gd as HTMLDivElement;
              }}
            />
          </div>
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      <section className={cardCls}>
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg text-black font-semibold">
            Parsimony{root ? ` — root=${root.replace(/^h_/, "")}` : ""} · layer={LAYER_LABEL[layer]}
          </h2>
          <button
            type="button"
            onClick={() =>
              downloadFigure(
                "par",
                `ecdll_job-${job || "unknown"}_root-${(root || "na").replace(/[^A-Za-z0-9_-]+/g, "_")}_rep-${
                  rep ?? "na"
                }_layer-${layer}_parsimony`
              )
            }
            className="px-2 py-1 rounded bg-black text-white"
          >
            Download PNG
          </button>
        </div>
        {parTrace ? (
          <div className={plotFrameCls}>
            <Plot
              data={[parTrace as Partial<PlotData>]}
              layout={baseLayout}
              config={{ displayModeBar: false, responsive: true }}
              style={{ width: "100%", height: "100%" }}
              onInitialized={(_, gd) => {
                graphRefs.current.par = gd as HTMLDivElement;
              }}
              onUpdate={(_, gd) => {
                graphRefs.current.par = gd as HTMLDivElement;
              }}
            />
          </div>
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      <section className={cardCls}>
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg text-black font-semibold">
            HSS{root ? ` — root=${root.replace(/^h_/, "")}` : ""} · layer={LAYER_LABEL[layer]}
          </h2>
          <button
            type="button"
            onClick={() =>
              downloadFigure(
                "hss",
                `ecdll_job-${job || "unknown"}_root-${(root || "na").replace(/[^A-Za-z0-9_-]+/g, "_")}_rep-${
                  rep ?? "na"
                }_layer-${layer}_hss`
              )
            }
            className="px-2 py-1 rounded bg-black text-white"
          >
            Download PNG
          </button>
        </div>
        {hssTrace ? (
          <div className={plotFrameCls}>
            <Plot
              data={[hssTrace as Partial<PlotData>]}
              layout={baseLayout}
              config={{ displayModeBar: false, responsive: true }}
              style={{ width: "100%", height: "100%" }}
              onInitialized={(_, gd) => {
                graphRefs.current.hss = gd as HTMLDivElement;
              }}
              onUpdate={(_, gd) => {
                graphRefs.current.hss = gd as HTMLDivElement;
              }}
            />
          </div>
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      {err ? <div className="p-3 border-2 border-black rounded text-red-600">{err}</div> : null}
    </div>
  );
}
