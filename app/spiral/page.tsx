// app/spiral/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import dynamic from "next/dynamic";
import type { PlotParams } from "react-plotly.js";
import type { Data, Layout, Config, ScatterData, PlotlyHTMLElement } from "plotly.js";
import { buildExportPrefix } from "@/lib/exportName";

const Plot = dynamic<PlotParams>(() => import("react-plotly.js"), { ssr: false });

type MethodName = "Dirichlet" | "Parsimony" | "SSH";

type SpiralRow = {
  root: string;
  rep: number;
  ll_init: number | null;
  ecd_ll_first: number | null;
  ecd_ll_final: number | null;
  ll_final: number | null;
};

type ApiSpiral = {
  job_id: string;
  method: MethodName;
  rows: SpiralRow[];
  reps: number[];
};

type ApiErr = { error: string };

// 17 internal nodes: h_21..h_37
const NODES = Array.from({ length: 17 }, (_, i) => `h_${21 + i}`);

// Sequential color scales
const BLUES_SCALE: [number, string][] = [
  [0.0,  "#f7fbff"],
  [0.125,"#deebf7"],
  [0.25, "#c6dbef"],
  [0.375,"#9ecae1"],
  [0.5,  "#6baed6"],
  [0.625,"#4292c6"],
  [0.75, "#2171b5"],
  [0.875,"#08519c"],
  [1.0,  "#08306b"],
];

const REDS_SCALE: [number, string][] = [
  [0.0,  "#fff5f0"],
  [0.125,"#fee0d2"],
  [0.25, "#fcbba1"],
  [0.375,"#fc9272"],
  [0.5,  "#fb6a4a"],
  [0.625,"#ef3b2c"],
  [0.75, "#cb181d"],
  [0.875,"#a50f15"],
  [1.0,  "#67000d"],
];

// ---------- small helpers ----------
function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isApiErr(x: unknown): x is ApiErr {
  return isRecord(x) && typeof x["error"] === "string";
}
function isApiSpiral(x: unknown): x is ApiSpiral {
  return (
    isRecord(x) &&
    typeof x["job_id"] === "string" &&
    (x["method"] === "Dirichlet" || x["method"] === "Parsimony" || x["method"] === "SSH") &&
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
function sampleReps(all: number[], k: number): number[] {
  if (all.length <= k) return [...all];
  const arr = [...all];
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [arr[i], arr[j]] = [arr[j], arr[i]];
  }
  return arr.slice(0, k).sort((a, b) => a - b);
}
const clamp01 = (x: number) => Math.max(0, Math.min(1, x));
const lerp = (a: number, b: number, t: number) => a + (b - a) * t;
function hexToRgb(hex: string) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex.trim());
  if (!m) return { r: 255, g: 255, b: 255 };
  const n = parseInt(m[1], 16);
  return { r: (n >> 16) & 255, g: (n >> 8) & 255, b: n & 255 };
}
function lerpRgb(a: {r:number;g:number;b:number}, b: {r:number;g:number;b:number}, t: number) {
  return {
    r: Math.round(lerp(a.r, b.r, t)),
    g: Math.round(lerp(a.g, b.g, t)),
    b: Math.round(lerp(a.b, b.b, t)),
  };
}
// colorScale as [[stop, "#hex"], ...], t in [0..1]
function lerpColors(scale: [number, string][], t: number) {
  const tt = clamp01(t);
  for (let i = 0; i < scale.length - 1; i++) {
    const [t0, c0] = scale[i]; const [t1, c1] = scale[i + 1];
    if (tt >= t0 && tt <= t1) {
      const local = (tt - t0) / Math.max(1e-9, (t1 - t0));
      const a = hexToRgb(c0), b = hexToRgb(c1);
      const m = lerpRgb(a, b, local);
      return `#${((1 << 24) + (m.r << 16) + (m.g << 8) + m.b).toString(16).slice(1)}`;
    }
  }
  return scale[tt < scale[0][0] ? 0 : scale.length - 1][1];
}

// ---------- node overlay (markers stay colored by metric, with colorbar) ----------
type NodeOverlayArgs = {
  coordFor: (node: string, rank: number) => { x: number; y: number };
  byRootRep: Map<string, Map<number, SpiralRow>>;
  selectedReps: number[];
  repRanks: Map<number, number>;
  colorScale: [number, string][];
  nodeColorMetric: (row: SpiralRow) => number | null;
  size?: number;
};
function buildNodeOverlayTrace(cfg: NodeOverlayArgs): Partial<ScatterData> {
  const { coordFor, byRootRep, selectedReps, repRanks, colorScale, nodeColorMetric, size = 10 } = cfg;

  // cmin/cmax
  const values: number[] = [];
  for (const node of NODES) {
    const perRep = byRootRep.get(node);
    if (!perRep) continue;
    for (const rep of selectedReps) {
      const row = perRep.get(rep);
      const v = row ? nodeColorMetric(row) : null;
      if (v != null && Number.isFinite(v)) values.push(v);
    }
  }
  let cMin = Math.min(...values), cMax = Math.max(...values);
  if (!values.length || !Number.isFinite(cMin) || !Number.isFinite(cMax) || cMin === cMax) { cMin = 0; cMax = 1; }

  const xs: number[] = [], ys: number[] = [], cs: number[] = [], texts: string[] = [];
  for (const node of NODES) {
    const perRep = byRootRep.get(node);
    if (!perRep) continue;
    for (const rep of selectedReps) {
      const row = perRep.get(rep);
      if (!row) continue;
      const rank = repRanks.get(rep) ?? 0;
      const { x, y } = coordFor(node, rank);
      const v = nodeColorMetric(row);
      if (v == null || !Number.isFinite(v)) continue;
      xs.push(x); ys.push(y); cs.push(v);
      texts.push(`${node}, rep ${rep}: ${v.toFixed(3)}`);
    }
  }
  return {
    type: "scatter",
    mode: "markers",
    x: xs, y: ys, text: texts, hoverinfo: "text",
    marker: {
      size,
      color: cs,
      colorscale: colorScale,
      cmin: cMin,
      cmax: cMax,
      reversescale: false,
      showscale: true,
      colorbar: {
        title: { text: "", font: { color: "white" } },
        tickfont: { color: "white" },
        thickness: 12,
        tickformat: ".3f",
      },
      line: { color: "white", width: 1 },
    },
    name: "nodes",
    showlegend: false,
  };
}

// ---------- small white center disk ----------
function buildCoreDot(radius: number, segments = 48): Partial<ScatterData> {
  const xs: number[] = [], ys: number[] = [];
  for (let i = 0; i <= segments; i++) {
    const th = (i / segments) * 2 * Math.PI;
    xs.push(radius * Math.cos(th));
    ys.push(radius * Math.sin(th));
  }
  return {
    type: "scatter",
    mode: "lines",
    x: xs, y: ys,
    line: { width: 0 },
    fill: "toself",
    fillcolor: "white",
    hoverinfo: "skip",
    showlegend: false,
    name: "core",
  };
}

// ---------- petal bands with inner↔outer gradient ----------
type PetalGradientArgs = {
  coordFor: (node: string, rank: number) => { x: number; y: number };
  byRootRep: Map<string, Map<number, SpiralRow>>;
  selectedReps: number[];
  repRanks: Map<number, number>;
  colorScale: [number, string][];
  nodeColorMetric: (row: SpiralRow) => number | null;

  petalWidthFrac?: number; // fraction of arm gap
  petalPower?: number;     // >1 => rounder tips
  samples?: number;        // angular samples per edge
  slicesPerBand?: number;  // radial slices to approximate gradient
  coreDotRadius?: number;  // white center disc radius
};
function buildPetalBandsGradient(cfg: PetalGradientArgs): Partial<ScatterData>[] {
  const {
    coordFor, byRootRep, selectedReps, repRanks,
    colorScale, nodeColorMetric,
    petalWidthFrac = 0.9, petalPower = 2.6,
    samples = 72, slicesPerBand = 4, coreDotRadius = 0.06,
  } = cfg;

  const traces: Partial<ScatterData>[] = [];
  traces.push(buildCoreDot(coreDotRadius, 48));
  if (!selectedReps.length) return traces;

  // global cmin/cmax
  const allVals: number[] = [];
  for (const node of NODES) {
    const perRep = byRootRep.get(node);
    if (!perRep) continue;
    for (const rep of selectedReps) {
      const row = perRep.get(rep);
      const v = row ? nodeColorMetric(row) : null;
      if (v != null && Number.isFinite(v)) allVals.push(v);
    }
  }
  let cMin = Math.min(...allVals), cMax = Math.max(...allVals);
  if (!allVals.length || !Number.isFinite(cMin) || !Number.isFinite(cMax) || cMin === cMax) { cMin = 0; cMax = 1; }

  // geometry
  const nArms = NODES.length;
  const dTheta = (2 * Math.PI) / nArms;
  const halfWidth = (dTheta * petalWidthFrac) / 2;

  // ranks for selected reps (stable ascending)
  const sortedReps = [...selectedReps].sort((a, b) => a - b);
  const m = sortedReps.length;

  // derive ring radii from coordFor (consistent with your spiral layout)
  function radiusFor(node: string, rank: number) {
    const p = coordFor(node, rank);
    return Math.hypot(p.x, p.y);
  }

  for (let i = 0; i < nArms; i++) {
    const node = NODES[i];
    const perRep = byRootRep.get(node);
    if (!perRep) continue;

    const theta0 = i * dTheta; // center angle for the petal

    const rowsByRank: (SpiralRow | null)[] = [];
    const rByRank: number[] = [];
    for (let rnk = 0; rnk < m; rnk++) {
      const rep = sortedReps[rnk];
      const row = perRep.get(rep) ?? null;
      rowsByRank.push(row);
      rByRank.push(radiusFor(node, repRanks.get(rep) ?? 0));
    }

    for (let k = 0; k < m - 1; k++) {
      const rowIn = rowsByRank[k];
      const rowOut = rowsByRank[k + 1];
      if (!rowIn && !rowOut) continue;

      const vIn = rowIn ? nodeColorMetric(rowIn) : null;
      const vOut = rowOut ? nodeColorMetric(rowOut) : null;

      const tIn = clamp01(((vIn ?? vOut ?? (cMin + cMax) / 2) - cMin) / Math.max(1e-9, (cMax - cMin)));
      const tOut = clamp01(((vOut ?? vIn ?? (cMin + cMax) / 2) - cMin) / Math.max(1e-9, (cMax - cMin)));
      const cIn = hexToRgb(lerpColors(colorScale, tIn));
      const cOut = hexToRgb(lerpColors(colorScale, tOut));

      const rIn = rByRank[k];
      const rOut = rByRank[k + 1];

      for (let s = 0; s < slicesPerBand; s++) {
        const a = s / slicesPerBand;
        const b = (s + 1) / slicesPerBand;

        const rLoBase = lerp(rIn, rOut, a);
        const rHiBase = lerp(rIn, rOut, b);

        const tBlend = (a + b) * 0.5;
        const cMid = lerpRgb(cIn, cOut, tBlend);
        const fill = `rgba(${cMid.r},${cMid.g},${cMid.b},0.95)`;

        const X: number[] = [];
        const Y: number[] = [];

        // Outer curve
        for (let j = 0; j <= samples; j++) {
          const u = (j / samples) * 2 - 1; // -1..+1
          const ang = -(theta0 + u * halfWidth);
          const sShape = Math.pow(Math.cos((Math.PI / 2) * u), petalPower);
          const r = coreDotRadius + (rHiBase - coreDotRadius) * sShape;
          X.push(r * Math.cos(ang));
          Y.push(r * Math.sin(ang));
        }
        // Inner curve (reverse)
        for (let j = samples; j >= 0; j--) {
          const u = (j / samples) * 2 - 1;
          const ang = -(theta0 + u * halfWidth);
          const sShape = Math.pow(Math.cos((Math.PI / 2) * u), petalPower);
          const r = coreDotRadius + (rLoBase - coreDotRadius) * sShape;
          X.push(r * Math.cos(ang));
          Y.push(r * Math.sin(ang));
        }

        traces.push({
          type: "scatter",
          mode: "lines",
          x: X, y: Y,
          hoverinfo: "skip",
          line: { width: 0 },        // seamless bands
          fill: "toself",
          fillcolor: fill,
          showlegend: false,
          name: `${node}-band-${k}-${s}`,
        });
      }
    }
  }

  return traces;
}

// ---------- petals + nodes wrapper ----------
type MetricCfg = {
  title: string;
  values: (r: SpiralRow) => number | null;
  colorScale: [number, string][];
};
type PetalOpts = {
  petalWidthFrac?: number;
  petalPower?: number;
  samples?: number;
  slicesPerBand?: number;
  coreDotRadius?: number;
};
function buildPetalBandsFromReps(
  metric: MetricCfg,
  coordFor: (node: string, rank: number) => { x: number; y: number },
  byRootRep: Map<string, Map<number, SpiralRow>>,
  selectedReps: number[],
  repRanks: Map<number, number>,
  opts?: PetalOpts
): Partial<ScatterData>[] {
  const petals = buildPetalBandsGradient({
    coordFor,
    byRootRep,
    selectedReps,
    repRanks,
    colorScale: metric.colorScale,
    nodeColorMetric: metric.values,
    petalWidthFrac: opts?.petalWidthFrac ?? 0.9,
    petalPower:     opts?.petalPower     ?? 2.6,
    samples:        opts?.samples        ?? 72,
    slicesPerBand:  opts?.slicesPerBand  ?? 4,
    coreDotRadius:  opts?.coreDotRadius  ?? 0.06,
  });

  const nodes = buildNodeOverlayTrace({
    coordFor,
    byRootRep,
    selectedReps,
    repRanks,
    colorScale: metric.colorScale,
    nodeColorMetric: metric.values,
    size: 10,
  });

  return [...petals, nodes];
}

// ---------- page ----------
export default function SpiralPage() {
  const [job, setJob] = useState<string>("");
  const [method, setMethod] = useState<MethodName>("Dirichlet");
  const [rows, setRows] = useState<SpiralRow[]>([]);
  const [repsAll, setRepsAll] = useState<number[]>([]);
  const [selectedReps, setSelectedReps] = useState<number[]>([]);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // read selected job (same key used across the app)
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
    setLoading(true); setErr(null);
    try {
      const u = new URL("/api/spiral", window.location.origin);
      u.searchParams.set("job", job);
      u.searchParams.set("method", method);
      const res = await fetch(u.toString(), { cache: "no-store" });

      const ct = res.headers.get("content-type") || "";
      if (!ct.includes("application/json")) {
        const body = await res.text().catch(() => "");
        throw new Error(
          `Expected JSON from ${u.pathname}, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`
        );
      }

      const json: unknown = await res.json();
      if (!res.ok || isApiErr(json)) throw new Error(isApiErr(json) ? json.error : `HTTP ${res.status}`);
      if (!isApiSpiral(json)) throw new Error("Invalid /api/spiral response");

      setRows(json.rows);
      setRepsAll(json.reps);
      setSelectedReps(sampleReps(json.reps, 5));
    } catch (e: unknown) {
      setErr(e instanceof Error ? e.message : "Failed to load");
      setRows([]); setRepsAll([]); setSelectedReps([]);
    } finally {
      setLoading(false);
    }
  }, [job, method]);

  useEffect(() => { if (job) void load(); }, [job, method, load]);

  const onResample = useCallback(() => {
    setSelectedReps(sampleReps(repsAll, 5));
  }, [repsAll]);

  // filter to selected reps
  const filtered = useMemo(() => {
    if (!rows.length || !selectedReps.length) return [] as SpiralRow[];
    const set = new Set(selectedReps);
    return rows.filter((r) => set.has(r.rep));
  }, [rows, selectedReps]);

  // root -> rep -> row
  const byRootRep = useMemo(() => {
    const map = new Map<string, Map<number, SpiralRow>>();
    for (const r of filtered) {
      if (!map.has(r.root)) map.set(r.root, new Map());
      map.get(r.root)!.set(r.rep, r);
    }
    return map;
  }, [filtered]);

  // ---------- spiral geometry for base radii (clockwise) ----------
  type XY = { x: number; y: number };
  type CoordKey = `${string}#${number}`;
  const coords = useMemo(() => {
    const out = new Map<CoordKey, XY>();
    const nArms = NODES.length;
    const dTheta = (2 * Math.PI) / nArms;
    const r0 = 0.6;          // base radius
    const dr = 0.35;         // per-rank radial step
    const twistDeg = -16;    // constant twist (deg/rank)
    const dTwist = (Math.PI / 180) * twistDeg;
    NODES.forEach((node, armIdx) => {
      const base = armIdx * dTheta;
      for (let k = 0; k < 5; k++) {
        const r = r0 + k * dr;
        const th = -(base + k * dTwist); // clockwise
        out.set(`${node}#${k}`, { x: r * Math.cos(th), y: r * Math.sin(th) });
      }
    });
    return out;
  }, []);

  const coordFor = useCallback(
    (node: string, rank: number): XY => coords.get(`${node}#${rank}`) ?? { x: 0, y: 0 },
    [coords]
  );

  // stable ranks 0..n-1 for selected reps
  const repRanks = useMemo(() => {
    const sorted = [...selectedReps].sort((a, b) => a - b);
    const r = new Map<number, number>();
    sorted.forEach((rep, idx) => r.set(rep, idx));
    return r;
  }, [selectedReps]);

  // ----- traces for each panel (petals + nodes) -----
  const initTraces = useMemo<Partial<ScatterData>[]>(() => {
    return buildPetalBandsFromReps(
      { title: "ll_init", values: (r) => r.ll_init, colorScale: BLUES_SCALE },
      coordFor, byRootRep, selectedReps, repRanks,
      { petalWidthFrac: 0.95, petalPower: 2.4, samples: 96, slicesPerBand: 5, coreDotRadius: 0.075 }
    );
  }, [coordFor, byRootRep, selectedReps, repRanks]);

  const finalTraces = useMemo<Partial<ScatterData>[]>(() => {
    return buildPetalBandsFromReps(
      { title: "ll_final", values: (r) => r.ll_final, colorScale: REDS_SCALE },
      coordFor, byRootRep, selectedReps, repRanks,
      { petalWidthFrac: 0.95, petalPower: 2.4, samples: 96, slicesPerBand: 5, coreDotRadius: 0.075 }
    );
  }, [coordFor, byRootRep, selectedReps, repRanks]);

  const ecdFirstTraces = useMemo<Partial<ScatterData>[]>(() => {
    return buildPetalBandsFromReps(
      { title: "ecd_ll_first", values: (r) => r.ecd_ll_first, colorScale: BLUES_SCALE },
      coordFor, byRootRep, selectedReps, repRanks,
      { petalWidthFrac: 0.95, petalPower: 2.4, samples: 96, slicesPerBand: 5, coreDotRadius: 0.075 }
    );
  }, [coordFor, byRootRep, selectedReps, repRanks]);

  const ecdFinalTraces = useMemo<Partial<ScatterData>[]>(() => {
    return buildPetalBandsFromReps(
      { title: "ecd_ll_final", values: (r) => r.ecd_ll_final, colorScale: REDS_SCALE },
      coordFor, byRootRep, selectedReps, repRanks,
      { petalWidthFrac: 0.95, petalPower: 2.4, samples: 96, slicesPerBand: 5, coreDotRadius: 0.075 }
    );
  }, [coordFor, byRootRep, selectedReps, repRanks]);

  // ----- export helpers -----
  const initGDRef = useRef<PlotlyHTMLElement | null>(null);
  const finalGDRef = useRef<PlotlyHTMLElement | null>(null);
  const ecdFirstGDRef = useRef<PlotlyHTMLElement | null>(null);
  const ecdFinalGDRef = useRef<PlotlyHTMLElement | null>(null);

  const exportPrefix = useMemo(
    () => buildExportPrefix({ job, twistDeg: -16 }),
    [job]
  );

  type ToImageOpts = { format?: "png" | "svg" | "jpeg" | "webp"; height?: number; width?: number; scale?: number };
  const toImage = useCallback(async (gd: PlotlyHTMLElement, opts: ToImageOpts) => {
    const w = window as unknown as { Plotly?: { toImage: (el: PlotlyHTMLElement, o: ToImageOpts) => Promise<string> } };
    if (!w.Plotly?.toImage) throw new Error("Plotly.toImage not available");
    return w.Plotly.toImage(gd, opts);
  }, []);
  const downloadURI = (uri: string, name: string) => {
    const a = document.createElement("a");
    a.href = uri; a.download = name; document.body.appendChild(a); a.click(); a.remove();
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

  // ---------- UI ----------
  return (
    <div className="min-h-screen bg-black text-white p-6 max-w-[1400px] mx-auto space-y-4">
      <h1 className="text-2xl font-bold">Flower Plots</h1>

      <div className="flex flex-wrap gap-3 items-end">
        <div className="text-sm text-gray-300">
          Job: <span className="font-mono">{job || "(none)"}</span>
        </div>

        {/* Method selector */}
        <div className="flex gap-2 ml-4">
          {(["Dirichlet", "Parsimony", "SSH"] as const).map((m) => (
            <button
              key={m}
              type="button"
              onClick={() => setMethod(m)}
              className={`px-3 py-2 rounded border ${
                method === m ? "bg-white text-black" : "bg-gray-800 text-white hover:bg-yellow-700"
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
            className="px-4 py-2 rounded bg-white text-black hover:bg-gray-100 disabled:opacity-50"
            disabled={!repsAll.length || loading}
          >
            Resample 5 reps
          </button>
          <button
            type="button"
            onClick={() => void downloadAllPlots()}
            className="px-4 py-2 rounded bg-white text-black hover:bg-gray-100"
            title="Download ll_init, ll_final, ecd_ll_first, ecd_ll_final as PNGs"
          >
            Download all 4
          </button>
        </div>
      </div>

      {loading ? (
        <div className="p-4 border border-gray-700 rounded bg-black text-white">Loading…</div>
      ) : err ? (
        <div className="p-4 border border-gray-700 rounded text-red-400 bg-black">{err}</div>
      ) : !rows.length || !selectedReps.length ? (
        <div className="p-4 border border-gray-700 rounded bg-black text-white">No data</div>
      ) : (
        <>
          {/* Row 1 */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <ChartCard title={`${method}: ll_init`} traces={initTraces} onReady={(gd) => (initGDRef.current = gd)} />
            <ChartCard title={`${method}: ll_final`} traces={finalTraces} onReady={(gd) => (finalGDRef.current = gd)} />
          </div>
          {/* Row 2 */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <ChartCard
              title={`${method}: ecd_ll_first`}
              traces={ecdFirstTraces}
              onReady={(gd) => (ecdFirstGDRef.current = gd)}
            />
            <ChartCard
              title={`${method}: ecd_ll_final`}
              traces={ecdFinalTraces}
              onReady={(gd) => (ecdFinalGDRef.current = gd)}
            />
          </div>
        </>
      )}
    </div>
  );
}

// ---------- card ----------
function ChartCard({
  title,
  traces,
  onReady,
  centerRadius = 0.08,
}: {
  title: string;
  traces: Partial<Data>[];
  onReady?: (el: PlotlyHTMLElement | null) => void;
  centerRadius?: number;
}) {
  const layout: Partial<Layout> = {
    title: { text: title, font: { color: "white" } },
    xaxis: { visible: false, scaleanchor: "y", color: "white" },
    yaxis: { visible: false, color: "white" },
    margin: { l: 20, r: 60, t: 40, b: 20 },
    height: 520,
    showlegend: false,
    paper_bgcolor: "black",
    plot_bgcolor: "black",
    shapes: [
      {
        type: "circle",
        xref: "x", yref: "y",
        x0: -centerRadius, y0: -centerRadius,
        x1:  centerRadius, y1:  centerRadius,
        line: { width: 0 },
        fillcolor: "white",
        layer: "above",
        opacity: 1,
      },
    ],
  };
  const config: Partial<Config> = { displayModeBar: false, responsive: true };

  return (
    <div className="bg-black border border-gray-700 rounded p-2">
      <Plot
        data={traces}
        layout={layout}
        config={config}
        style={{ width: "100%", height: "100%" }}
        onInitialized={(_figure: unknown, gd: PlotlyHTMLElement) => onReady?.(gd)}
        onUpdate={(_figure: unknown, gd: PlotlyHTMLElement) => onReady?.(gd)}
      />
    </div>
  );
}
