// app/dosa/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import dynamic from "next/dynamic";
import type { PlotParams } from "react-plotly.js";
import type { Data, Layout, Config, ScatterData, PlotlyHTMLElement } from "plotly.js";
import { buildExportPrefix } from "@/lib/exportName";

const Plot = dynamic<PlotParams>(() => import("react-plotly.js"), { ssr: false });

type MethodName = "Dirichlet" | "Parsimony" | "SSH";

type RhodoneaRow = {
  root: string;
  rep: number;
  ll_init: number | null;
  ecd_ll_first: number | null;
  ecd_ll_final: number | null;
  ll_final: number | null;
};

type ApiRhodonea = {
  job_id: string;
  method: string;
  rows: RhodoneaRow[];
  reps: number[];
};

type ApiErr = { error: string };

// 17 internal nodes: h_21..h_37
const NODES = Array.from({ length: 17 }, (_, i) => `h_${21 + i}`);

// Sequential color scales
const BLUES_SCALE: [number, string][] = [
  [0.0,   "#6baed6"],
  [0.125, "#57a0ce"],
  [0.25,  "#4292c6"],
  [0.375, "#3282be"],
  [0.5,   "#2171b5"],
  [0.625, "#1461a8"],
  [0.75,  "#08519c"],
  [0.875, "#084084"],
  [1.0,   "#08306b"],
];

const REDS_SCALE: [number, string][] = [
  [0.0,   "#fb6a4a"],
  [0.125, "#f5523b"],
  [0.25,  "#ef3b2c"],
  [0.375, "#dd2a24"],
  [0.5,   "#cb181d"],
  [0.625, "#b81419"],
  [0.75,  "#a50f15"],
  [0.875, "#860811"],
  [1.0,   "#67000d"],
];

/* -------------------- Utils -------------------- */
function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isApiErr(x: unknown): x is ApiErr {
  return isRecord(x) && typeof x["error"] === "string";
}
function isApiRhodonea(x: unknown): x is ApiRhodonea {
  return (
    isRecord(x) &&
    typeof x["job_id"] === "string" &&
    typeof x["method"] === "string" &&
    Array.isArray(x["rows"]) &&
    Array.isArray(x["reps"])
  );
}
function extent(nums: number[]): [number, number] {
  if (!nums.length) return [0, 1];
  let mn = nums[0], mx = nums[0];
  for (const v of nums) { if (v < mn) mn = v; if (v > mx) mx = v; }
  if (mn === mx) return [mn - 1, mx + 1];
  return [mn, mx];
}
function clamp01(x: number) { return Math.max(0, Math.min(1, x)); }
function colorFromScale(
  v: number,
  vmin: number,
  vmax: number,
  scale: [number, string][]
): string {
  if (!Number.isFinite(v)) return "rgb(255,255,255)";
  const t = vmax === vmin ? 0.5 : clamp01((v - vmin) / (vmax - vmin));
  for (let i = 1; i < scale.length; i++) {
    const [t1, c1] = scale[i - 1];
    const [t2, c2] = scale[i];
    if (t <= t2) {
      const u = clamp01((t - t1) / Math.max(1e-9, t2 - t1));
      const parse = (s: string) => {
        const m = s.match(/#([0-9a-f]{6})/i);
        if (!m) return [255, 255, 255] as const;
        const n = parseInt(m[1], 16);
        return [(n >> 16) & 255, (n >> 8) & 255, n & 255] as const;
      };
      const [r1, g1, b1] = parse(c1);
      const [r2, g2, b2] = parse(c2);
      const r = Math.round(r1 + u * (r2 - r1));
      const g = Math.round(g1 + u * (g2 - g1)); // ✅ correct channel
      const b = Math.round(b1 + u * (b2 - b1));
      return `rgb(${r},${g},${b})`;
    }
  }
  return scale[scale.length - 1][1];
}

/* -------------------- Page -------------------- */
export default function RhodoneaPage() {
  const [job, setJob] = useState<string>("");
  const [method, setMethod] = useState<MethodName>("Dirichlet");
  const [rows, setRows] = useState<RhodoneaRow[]>([]);
  const [repsAll, setRepsAll] = useState<number[]>([]);
  const [selectedReps, setSelectedReps] = useState<number[]>([]);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // geometry constants for dosa layout
  const r0 = 0.6;           // inner ring radius
  const dr = 0.35;          // ring spacing
  const coreRadius = 0.085; // center disk & dosa base
  const turns = 1.25;       // total windings from base to tip
  const ribbonHalfWidth = 0.035; // angular half-thickness of ribbon (radians)

  // read selected job
  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId");
      if (saved) setJob(saved);
    } catch {/* ignore */}
  }, []);

  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Choose a job on Precomputed Results first.");
      setRows([]); setRepsAll([]); setSelectedReps([]);
      return;
    }
    setLoading(true);
    setErr(null);
    try {
      const u = new URL("/api/dosa", window.location.origin);
      u.searchParams.set("job", job);
      u.searchParams.set("method", method);
      const res = await fetch(u.toString(), { cache: "no-store" });

      const raw = await res.text();
      let json: unknown;
      try {
        json = raw ? JSON.parse(raw) : {};
      } catch {
        throw new Error(
          `Expected JSON from ${u.pathname}, got invalid body (HTTP ${res.status}). ${raw.slice(0, 200)}`
        );
      }
      if (!res.ok || isApiErr(json)) {
        throw new Error(isApiErr(json) ? json.error : `HTTP ${res.status}`);
      }
      if (!isApiRhodonea(json)) {
        throw new Error("Invalid /api/dosa response");
      }

      setRows(json.rows);
      setRepsAll(json.reps);
      setSelectedReps(sampleReps(json.reps, Math.min(5, json.reps.length)));
    } catch (e: unknown) {
      setErr(e instanceof Error ? e.message : "Failed to load");
      setRows([]); setRepsAll([]); setSelectedReps([]);
    } finally {
      setLoading(false);
    }
  }, [job, method]);

  useEffect(() => { if (job) void load(); }, [job, method, load]);

  const onResample = useCallback(() => {
    setSelectedReps(sampleReps(repsAll, Math.min(5, repsAll.length)));
  }, [repsAll]);

  // ----- Build rep ranking across ALL rows (for outermost detection) -----
  const allRepsSorted = useMemo(
    () => Array.from(new Set(rows.map((r) => r.rep))).sort((a, b) => a - b),
    [rows]
  );
  const repRankAll = useMemo(() => {
    const m = new Map<number, number>();
    allRepsSorted.forEach((rep, idx) => m.set(rep, idx));
    return m;
  }, [allRepsSorted]);

  // Map: root -> (rep -> row) using ALL rows
  const byRootRepAll = useMemo(() => {
    const map = new Map<string, Map<number, RhodoneaRow>>();
    for (const r of rows) {
      if (!map.has(r.root)) map.set(r.root, new Map());
      map.get(r.root)!.set(r.rep, r);
    }
    return map;
  }, [rows]);

  // ---------- Label trace (outermost) ----------
  const labelTrace = useMemo(() => {
    return buildOuterLabelsTrace(byRootRepAll, repRankAll, {
      r0, dr,
      labelOffset: 0.06, // radial push
      tOffset: 0.03,     // tangential nudge
    });
  }, [byRootRepAll, repRankAll, r0, dr]);

  // ---------- Traces (Dosa ribbons + selected markers + labels) ----------
  const initTraces = useMemo<Partial<ScatterData>[]>(() => {
    if (!rows.length || !selectedReps.length) return [];
    const avgMap = makeNodeAvgMap(byRootRepAll, selectedReps, (r) => r.ll_init);
    const avgs = Array.from(avgMap.values());
    const [mn, mx] = extent(avgs);
    const dosas = buildDosaRibbonTraces(
      { title: "ll_init", colorScale: BLUES_SCALE, cmin: mn, cmax: mx },
      byRootRepAll, repRankAll, avgMap,
      { turns, ribbonHalfWidth, samples: 220, coreRadius, r0, dr }
    );
    const markers = buildMarkersTraceDosa(
      { title: "ll_init", values: (r) => r.ll_init, colorScale: BLUES_SCALE, cmin: mn, cmax: mx },
      byRootRepAll, repRankAll, selectedReps,
      { coreRadius, r0, dr, turns }
    );
    return [...dosas, markers, ...(labelTrace ? [labelTrace as Partial<ScatterData>] : [])];
  }, [rows, byRootRepAll, repRankAll, selectedReps, coreRadius, r0, dr, turns, ribbonHalfWidth, labelTrace]);

  const finalTraces = useMemo<Partial<ScatterData>[]>(() => {
    if (!rows.length || !selectedReps.length) return [];
    const avgMap = makeNodeAvgMap(byRootRepAll, selectedReps, (r) => r.ll_final);
    const avgs = Array.from(avgMap.values());
    const [mn, mx] = extent(avgs);
    const dosas = buildDosaRibbonTraces(
      { title: "ll_final", colorScale: REDS_SCALE, cmin: mn, cmax: mx },
      byRootRepAll, repRankAll, avgMap,
      { turns, ribbonHalfWidth, samples: 220, coreRadius, r0, dr }
    );
    const markers = buildMarkersTraceDosa(
      { title: "ll_final", values: (r) => r.ll_final, colorScale: REDS_SCALE, cmin: mn, cmax: mx },
      byRootRepAll, repRankAll, selectedReps,
      { coreRadius, r0, dr, turns }
    );
    return [...dosas, markers, ...(labelTrace ? [labelTrace as Partial<ScatterData>] : [])];
  }, [rows, byRootRepAll, repRankAll, selectedReps, coreRadius, r0, dr, turns, ribbonHalfWidth, labelTrace]);

  const ecdFirstTraces = useMemo<Partial<ScatterData>[]>(() => {
    if (!rows.length || !selectedReps.length) return [];
    const avgMap = makeNodeAvgMap(byRootRepAll, selectedReps, (r) => r.ecd_ll_first);
    const avgs = Array.from(avgMap.values());
    const [mn, mx] = extent(avgs);
    const dosas = buildDosaRibbonTraces(
      { title: "ecd_ll_first", colorScale: BLUES_SCALE, cmin: mn, cmax: mx },
      byRootRepAll, repRankAll, avgMap,
      { turns, ribbonHalfWidth, samples: 220, coreRadius, r0, dr }
    );
    const markers = buildMarkersTraceDosa(
      { title: "ecd_ll_first", values: (r) => r.ecd_ll_first, colorScale: BLUES_SCALE, cmin: mn, cmax: mx },
      byRootRepAll, repRankAll, selectedReps,
      { coreRadius, r0, dr, turns }
    );
    return [...dosas, markers, ...(labelTrace ? [labelTrace as Partial<ScatterData>] : [])];
  }, [rows, byRootRepAll, repRankAll, selectedReps, coreRadius, r0, dr, turns, ribbonHalfWidth, labelTrace]);

  const ecdFinalTraces = useMemo<Partial<ScatterData>[]>(() => {
    if (!rows.length || !selectedReps.length) return [];
    const avgMap = makeNodeAvgMap(byRootRepAll, selectedReps, (r) => r.ecd_ll_final);
    const avgs = Array.from(avgMap.values());
    const [mn, mx] = extent(avgs);
    const dosas = buildDosaRibbonTraces(
      { title: "ecd_ll_final", colorScale: REDS_SCALE, cmin: mn, cmax: mx },
      byRootRepAll, repRankAll, avgMap,
      { turns, ribbonHalfWidth, samples: 220, coreRadius, r0, dr }
    );
    const markers = buildMarkersTraceDosa(
      { title: "ecd_ll_final", values: (r) => r.ecd_ll_final, colorScale: REDS_SCALE, cmin: mn, cmax: mx },
      byRootRepAll, repRankAll, selectedReps,
      { coreRadius, r0, dr, turns }
    );
    return [...dosas, markers, ...(labelTrace ? [labelTrace as Partial<ScatterData>] : [])];
  }, [rows, byRootRepAll, repRankAll, selectedReps, coreRadius, r0, dr, turns, ribbonHalfWidth, labelTrace]);

  // ---------- Export helpers ----------
  const initGDRef = useRef<PlotlyHTMLElement | null>(null);
  const finalGDRef = useRef<PlotlyHTMLElement | null>(null);
  const ecdFirstGDRef = useRef<PlotlyHTMLElement | null>(null);
  const ecdFinalGDRef = useRef<PlotlyHTMLElement | null>(null);

  const exportPrefix = useMemo(() => buildExportPrefix({ job }), [job]);

  type ToImageOpts = { format?: "png" | "svg" | "jpeg" | "webp"; height?: number; width?: number; scale?: number };
  const toImage = useCallback(async (gd: PlotlyHTMLElement, opts: ToImageOpts) => {
    const w = window as unknown as { Plotly?: { toImage: (el: PlotlyHTMLElement, o: ToImageOpts) => Promise<string> } };
    if (!w.Plotly?.toImage) throw new Error("Plotly.toImage not available");
    return w.Plotly.toImage(gd, opts);
  }, []);
  const downloadURI = (uri: string, name: string) => {
    const a = document.createElement("a");
    a.href = uri; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
  };
  const downloadAllPlots = useCallback(async () => {
    const stamp = new Date().toISOString().replace(/[:.]/g, "").replace("T", "_").slice(0, 15);
    const grab = async (gd: PlotlyHTMLElement | null, metric: string) => {
      if (!gd) return;
      const url = await toImage(gd, { format: "png", scale: 2 });
      downloadURI(url, `${exportPrefix}__${method}__${metric}__t${stamp}.png`);
    };
    await Promise.all([
      grab(initGDRef.current, "ll_init"),
      grab(finalGDRef.current, "ll_final"),
      grab(ecdFirstGDRef.current, "ecd_ll_first"),
      grab(ecdFinalGDRef.current, "ecd_ll_final"),
    ]);
  }, [toImage, exportPrefix, method]);

  return (
    <div className="min-h-screen bg-white text-black p-6 max-w-[1400px] mx-auto space-y-4">
      <h1 className="text-2xl font-bold">Dosa Plots</h1>

      <div className="flex flex-wrap gap-3 items-end">
        <div className="text-sm text-gray-600">
          Job: <span className="font-mono">{job || "(none)"}</span>
        </div>

        <div className="flex gap-2 ml-4">
          {(["Dirichlet", "Parsimony", "SSH"] as const).map((m) => (
            <button
              key={m}
              type="button"
              onClick={() => setMethod(m)}
              className={`px-3 py-2 rounded border ${
                method === m
                  ? "bg-black text-white border-black"
                  : "bg-gray-100 text-black hover:bg-gray-200 border-gray-300"
              }`}
            >
              {m}
            </button>
          ))}
        </div>

        <div className="ml-auto flex gap-2">
          <button
            type="button"
            onClick={onResample}
            className="px-4 py-2 rounded bg-black text-white hover:opacity-90 disabled:opacity-50"
            disabled={!repsAll.length || loading}
          >
            Resample 5 reps
          </button>
          <button
            type="button"
            onClick={() => void downloadAllPlots()}
            className="px-4 py-2 rounded bg-black text-white hover:opacity-90"
            title="Download ll_init, ll_final, ecd_ll_first, ecd_ll_final as PNGs"
          >
            Download all 4
          </button>
        </div>
      </div>

      {loading ? (
        <div className="p-4 border border-gray-200 rounded bg-white text-black">Loading…</div>
      ) : err ? (
        <div className="p-4 border border-gray-200 rounded text-red-600 bg-white">{err}</div>
      ) : !rows.length || !selectedReps.length ? (
        <div className="p-4 border border-gray-200 rounded bg-white text-black">No data</div>
      ) : (
        <>
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <ChartCard title={`${method}: ll_init`} traces={initTraces} onReady={(gd) => (initGDRef.current = gd)} centerRadius={coreRadius} />
            <ChartCard title={`${method}: ll_final`} traces={finalTraces} onReady={(gd) => (finalGDRef.current = gd)} centerRadius={coreRadius} />
          </div>
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <ChartCard title={`${method}: ecd_ll_first`} traces={ecdFirstTraces} onReady={(gd) => (ecdFirstGDRef.current = gd)} centerRadius={coreRadius} />
            <ChartCard title={`${method}: ecd_ll_final`} traces={ecdFinalTraces} onReady={(gd) => (ecdFinalGDRef.current = gd)} centerRadius={coreRadius} />
          </div>
        </>
      )}
    </div>
  );

  function sampleReps(all: number[], k: number): number[] {
    if (k <= 0) return [];
    if (all.length <= k) return [...all].sort((a, b) => a - b);
    const arr = [...all];
    for (let i = arr.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [arr[i], arr[j]] = [arr[j], arr[i]];
    }
    return arr.slice(0, k).sort((a, b) => a - b);
  }
}

/* -------------------- Drawing -------------------- */

/** Average a metric over the SELECTED reps for each node. */
function makeNodeAvgMap(
  byRootRepAll: Map<string, Map<number, RhodoneaRow>>,
  selectedReps: number[],
  values: (r: RhodoneaRow) => number | null
): Map<string, number> {
  const m = new Map<string, number>();
  for (const node of NODES) {
    const perRep = byRootRepAll.get(node);
    if (!perRep) continue;
    let sum = 0, cnt = 0;
    for (const rep of selectedReps) {
      const row = perRep.get(rep);
      if (!row) continue;
      const v = values(row);
      if (v != null && Number.isFinite(v)) { sum += v; cnt++; }
    }
    if (cnt > 0) m.set(node, sum / cnt);
  }
  return m;
}

/** Build dosa ribbons (Archimedean) for each node arm. */
function buildDosaRibbonTraces(
  cfg: { title: string; colorScale: [number, string][]; cmin: number; cmax: number; },
  byRootRepAll: Map<string, Map<number, RhodoneaRow>>,
  repRankAll: Map<number, number>,
  avgMap: Map<string, number>,
  opts?: {
    turns?: number; ribbonHalfWidth?: number; samples?: number;
    coreRadius?: number; r0?: number; dr?: number;
  }
): Partial<ScatterData>[] {
  const traces: Partial<ScatterData>[] = [];
  const turns = Math.max(0, opts?.turns ?? 1.25);
  const ribbonHalfWidth = Math.max(0.001, opts?.ribbonHalfWidth ?? 0.035);
  const samples = Math.max(72, Math.min(1024, opts?.samples ?? 220));
  const coreRadius = Math.max(0.0, opts?.coreRadius ?? 0.085);
  const r0 = opts?.r0 ?? 0.6;
  const dr = opts?.dr ?? 0.35;

  const n = NODES.length;
  if (!n) return traces;

  const dThetaArm = (2 * Math.PI) / n;

  for (let i = 0; i < n; i++) {
    const node = NODES[i];
    const perRep = byRootRepAll.get(node);
    if (!perRep) continue;

    let outerRank = -1;
    let outerRep: number | null = null;
    for (const rep of perRep.keys()) {
      const rnk = repRankAll.get(rep);
      if (rnk == null) continue;
      if (rnk > outerRank) { outerRank = rnk; outerRep = rep; }
    }
    if (outerRank < 0 || outerRep == null) continue;

    const avg = avgMap.get(node);
    if (avg == null || !Number.isFinite(avg)) continue;

    const theta0 = -(i * dThetaArm);
    const Rmax = r0 + outerRank * dr;
    const span = Math.max(0, Rmax - coreRadius);
    if (span <= 0) continue;

    // Archimedean dosa: r(t) = core + span * t; theta(t) = theta0 + 2π * turns * t
    const leftX: number[] = [];
    const leftY: number[] = [];
    const rightX: number[] = [];
    const rightY: number[] = [];

    for (let j = 0; j <= samples; j++) {
      const t = j / samples;
      const r = coreRadius + span * t;
      const th = theta0 + (2 * Math.PI * turns) * t;

      // thin ribbon by offsetting ±ribbonHalfWidth in angle
      leftX.push(r * Math.cos(th - ribbonHalfWidth));
      leftY.push(r * Math.sin(th - ribbonHalfWidth));
      rightX.push(r * Math.cos(th + ribbonHalfWidth));
      rightY.push(r * Math.sin(th + ribbonHalfWidth));
    }

    // close polygon: left forward, right backward
    const xs = leftX.concat(rightX.slice().reverse());
    const ys = leftY.concat(rightY.slice().reverse());

    traces.push({
      type: "scatter",
      mode: "lines",
      x: xs,
      y: ys,
      fill: "toself",
      fillcolor: colorFromScale(avg, cfg.cmin, cfg.cmax, cfg.colorScale),
      // ⬇ better contrast on light bg
      line: { width: 0.6, color: "rgba(0,0,0,0.25)" },
      hovertemplate: `${node}, outermost rep ${outerRep}<br>${cfg.title} (avg of selected): ${avg.toFixed(3)}<extra></extra>`,
      showlegend: false,
      name: `${node}-dosa`,
    });
  }

  return traces;
}

/** Overlay replicate markers along the same dosa, spaced between base and tip. */
function buildMarkersTraceDosa(
  cfg: {
    title: string;
    values: (r: RhodoneaRow) => number | null;
    colorScale: [number, string][];
    cmin: number;
    cmax: number;
  },
  byRootRepAll: Map<string, Map<number, RhodoneaRow>>,
  repRankAll: Map<number, number>,
  selectedReps: number[],
  geometry: { coreRadius: number; r0: number; dr: number; turns?: number }
): Partial<ScatterData> {
  const xs: number[] = [];
  const ys: number[] = [];
  const colors: number[] = [];
  const texts: string[] = [];

  const nArms = NODES.length;
  const dThetaArm = (2 * Math.PI) / nArms;
  const turns = geometry.turns ?? 1.25;

  for (let i = 0; i < nArms; i++) {
    const node = NODES[i];
    const perRep = byRootRepAll.get(node);
    if (!perRep) continue;

    const repsHere = selectedReps.filter((rep) => perRep.has(rep));
    if (!repsHere.length) continue;

    let outerRank = -1;
    for (const rep of perRep.keys()) {
      const rnk = repRankAll.get(rep);
      if (rnk == null) continue;
      if (rnk > outerRank) outerRank = rnk;
    }
    if (outerRank < 0) continue;

    const theta0 = -(i * dThetaArm);
    const Rmax = geometry.r0 + outerRank * geometry.dr;
    const base = geometry.coreRadius + 0.02;
    const span = Math.max(0, Rmax - base);
    if (span <= 0) continue;

    repsHere.sort((a, b) => a - b);
    const m = repsHere.length;

    for (let k = 0; k < m; k++) {
      const frac = (k + 1) / (m + 1);      // (0,1) spacing
      const r = base + frac * span;        // radial
      const t = (r - base) / span;         // 0..1 param
      const th = theta0 + (2 * Math.PI * turns) * t;

      const row = perRep.get(repsHere[k])!;
      const v = cfg.values(row);
      if (v == null || !Number.isFinite(v)) continue;

      xs.push(r * Math.cos(th));
      ys.push(r * Math.sin(th));
      colors.push(v);
      texts.push(`${node}, rep ${repsHere[k]}<br>${cfg.title}: ${v.toFixed(3)}`);
    }
  }

  return {
    type: "scatter",
    mode: "markers",
    x: xs,
    y: ys,
    text: texts,
    hoverinfo: "text",
    marker: {
      size: 10,
      color: colors,
      colorscale: cfg.colorScale,
      reversescale: false,
      cmin: cfg.cmin,
      cmax: cfg.cmax,
      showscale: true,
      colorbar: {
        title: { text: cfg.title, font: { color: "black" } },
        tickfont: { color: "black" },
        thickness: 12,
        tickformat: ".3f",
      },
      line: { color: "black", width: 1 },
    },
    name: cfg.title,
    showlegend: false,
  };
}

/** Text labels placed just beyond each arm's tip, oriented radially. */
function buildOuterLabelsTrace(
  byRoot: Map<string, Map<number, RhodoneaRow>>,
  repRanks: Map<number, number>,
  geometry: { r0: number; dr: number; labelOffset?: number; tOffset?: number }
): Partial<ScatterData> | null {
  const n = NODES.length;
  if (!n) return null;

  const { r0, dr, labelOffset = 0.06, tOffset = 0.0 } = geometry;
  const dTheta = (2 * Math.PI) / n;

  const xs: number[] = [];
  const ys: number[] = [];
  const texts: string[] = [];

  for (let i = 0; i < n; i++) {
    const node = NODES[i];
    const perRep = byRoot.get(node);
    if (!perRep) continue;

    // outermost ring index (same as dosa tip)
    let outerRank = -1;
    for (const rep of perRep.keys()) {
      const rnk = repRanks.get(rep);
      if (rnk == null) continue;
      if (rnk > outerRank) outerRank = rnk;
    }
    if (outerRank < 0) continue;

    const theta0 = -(i * dTheta);
    const Rmax   = r0 + outerRank * dr;
    const Rlabel = Rmax + labelOffset;

    // base: radial position just beyond the tip
    let x = Rlabel * Math.cos(theta0);
    let y = Rlabel * Math.sin(theta0);

    // small tangential nudge; t̂ = (-sin θ, cos θ)
    if (tOffset !== 0) {
      x += tOffset * (-Math.sin(theta0));
      y += tOffset * ( Math.cos(theta0));
    }

    xs.push(x);
    ys.push(y);
    texts.push(node);
  }

  if (!xs.length) return null;

  return {
    type: "scatter",
    mode: "text",
    x: xs,
    y: ys,
    text: texts,
    textposition: "middle center",
    textfont: { color: "black", size: 11 },
    hoverinfo: "skip",
    showlegend: false,
    cliponaxis: false,
    name: "node labels",
  };
}

/** Plot card with a small center disk to anchor dosas visually. */
function ChartCard({
  title,
  traces,
  onReady,
  centerRadius = 0.085,
}: {
  title: string;
  traces: Partial<Data>[];
  onReady?: (el: PlotlyHTMLElement | null) => void;
  centerRadius?: number;
}) {
  const layout: Partial<Layout> = {
    title: { text: title, font: { color: "black" } },
    xaxis: { visible: false, scaleanchor: "y", color: "black" },
    yaxis: { visible: false, color: "black" },
    margin: { l: 20, r: 60, t: 40, b: 20 },
    height: 520,
    showlegend: false,
    paper_bgcolor: "white",
    plot_bgcolor: "white",
    shapes: [
      {
        type: "circle",
        xref: "x", yref: "y",
        x0: -centerRadius, y0: -centerRadius,
        x1:  centerRadius, y1:  centerRadius,
        line: { width: 0 },
        fillcolor: "black",
        layer: "below",   // ⬅ so it doesn't cover traces
        opacity: 1,
      },
    ],
  };
  const config: Partial<Config> = { displayModeBar: false, responsive: true };

  return (
    <div className="bg-white border border-gray-200 rounded p-2 shadow-sm">
      <Plot
        data={traces as Data[]}
        layout={layout}
        config={config}
        style={{ width: "100%", height: "100%" }}
        onInitialized={(_figure: unknown, gd: PlotlyHTMLElement) => onReady?.(gd)}
        onUpdate={(_figure: unknown, gd: PlotlyHTMLElement) => onReady?.(gd)}
      />
    </div>
  );
}
