// app/spiral/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import dynamic from "next/dynamic";
import type { PlotParams } from "react-plotly.js";
import type { Data, Layout, Config, ScatterData } from "plotly.js";

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

const NODES = Array.from({ length: 17 }, (_, i) => `h_${21 + i}`);

// Explicit light→dark sequential scales to guarantee orientation
const BLUES_SCALE: [number, string][] = [
  [0, "#f7fbff"],
  [0.125, "#deebf7"],
  [0.25, "#c6dbef"],
  [0.375, "#9ecae1"],
  [0.5, "#6baed6"],
  [0.625, "#4292c6"],
  [0.75, "#2171b5"],
  [0.875, "#08519c"],
  [1, "#08306b"],
];

const REDS_SCALE: [number, string][] = [
  [0, "#fff5f0"],
  [0.125, "#fee0d2"],
  [0.25, "#fcbba1"],
  [0.375, "#fc9272"],
  [0.5, "#fb6a4a"],
  [0.625, "#ef3b2c"],
  [0.75, "#cb181d"],
  [0.875, "#a50f15"],
  [1, "#67000d"],
];

// ------- type guards -------
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

export default function SpiralPage() {
  const [job, setJob] = useState<string>("");
  const [method, setMethod] = useState<MethodName>("Dirichlet");
  const [rows, setRows] = useState<SpiralRow[]>([]);
  const [repsAll, setRepsAll] = useState<number[]>([]);
  const [selectedReps, setSelectedReps] = useState<number[]>([]);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // User-controlled twist in degrees per rank (0 => radial)
  const [twistDeg, setTwistDeg] = useState<number>(18); // ≈ π/10 per-rank

  // read selected job (same key used across the app)
  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId");
      if (saved) setJob(saved);
    } catch {
      /* ignore */
    }
  }, []);

  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Choose a job on Precomputed Results first.");
      setRows([]);
      setRepsAll([]);
      setSelectedReps([]);
      return;
    }
    setLoading(true);
    setErr(null);
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
      setRows([]);
      setRepsAll([]);
      setSelectedReps([]);
    } finally {
      setLoading(false);
    }
  }, [job, method]);

  useEffect(() => {
    if (job) void load();
  }, [job, method, load]);

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

  // ---------- CLOCKWISE spiral geometry ----------
  type XY = { x: number; y: number };
  type CoordKey = `${string}#${number}`; // "h_21#rank"
  const coords = useMemo(() => {
    const out = new Map<CoordKey, XY>();
    const nArms = NODES.length; // 17
    const dTheta = (2 * Math.PI) / nArms; // arm spacing
    const r0 = 0.6;  // base radius
    const dr = 0.35; // per-rank radial step

    // degrees -> radians per rank
    const dTwist = (Math.PI / 180) * twistDeg;

    NODES.forEach((node, armIdx) => {
      const base = armIdx * dTheta;
      for (let k = 0; k < 5; k++) {
        const r = r0 + k * dr;
        const th = -(base + k * dTwist); // NEGATED => clockwise; 0 => radial
        out.set(`${node}#${k}`, { x: r * Math.cos(th), y: r * Math.sin(th) });
      }
    });
    return out;
  }, [twistDeg]);

  // stable rank 0..4 for whichever reps were selected
  const repRanks = useMemo(() => {
    const sorted = [...selectedReps].sort((a, b) => a - b);
    const r = new Map<number, number>();
    sorted.forEach((rep, idx) => r.set(rep, idx));
    return r;
  }, [selectedReps]);

  // ---- stable helpers ----
  const coordFor = useCallback(
    (node: string, rank: number): XY => {
      const key: CoordKey = `${node}#${rank}`;
      return coords.get(key) ?? { x: 0, y: 0 };
    },
    [coords]
  );

  const collectMetric = useCallback(
    (values: (row: SpiralRow) => number | null): number[] => {
      const out: number[] = [];
      for (const node of NODES) {
        const perRep = byRootRep.get(node);
        if (!perRep) continue;
        for (const rep of selectedReps) {
          const row = perRep.get(rep);
          if (!row) continue;
          const v = values(row);
          if (v != null && Number.isFinite(v)) out.push(v);
        }
      }
      return out;
    },
    [byRootRep, selectedReps]
  );

  // ----- Build each plot’s traces (dotted replicate links + markers with colorbar) -----
  const initTraces = useMemo<Partial<ScatterData>[]>(() => {
    const vals = collectMetric((r) => r.ll_init);
    const [mn, mx] = extent(vals);
    return [
      buildConnectionsTrace(coordFor, byRootRep, selectedReps, repRanks),
      buildMarkersTrace(
        {
          title: "ll_init",
          values: (r) => r.ll_init,
          colorScale: BLUES_SCALE, // light -> dark
          cmin: mn,
          cmax: mx,
        },
        coordFor,
        byRootRep,
        selectedReps,
        repRanks
      ),
    ];
  }, [collectMetric, coordFor, byRootRep, selectedReps, repRanks]);

  const finalTraces = useMemo<Partial<ScatterData>[]>(() => {
    const vals = collectMetric((r) => r.ll_final);
    const [mn, mx] = extent(vals);
    return [
      buildConnectionsTrace(coordFor, byRootRep, selectedReps, repRanks),
      buildMarkersTrace(
        {
          title: "ll_final",
          values: (r) => r.ll_final,
          colorScale: REDS_SCALE, // light -> dark
          cmin: mn,
          cmax: mx,
        },
        coordFor,
        byRootRep,
        selectedReps,
        repRanks
      ),
    ];
  }, [collectMetric, coordFor, byRootRep, selectedReps, repRanks]);

  const ecdFirstTraces = useMemo<Partial<ScatterData>[]>(() => {
    const vals = collectMetric((r) => r.ecd_ll_first);
    const [mn, mx] = extent(vals);
    return [
      buildConnectionsTrace(coordFor, byRootRep, selectedReps, repRanks),
      buildMarkersTrace(
        {
          title: "ecd_ll_first",
          values: (r) => r.ecd_ll_first,
          colorScale: BLUES_SCALE, // light -> dark
          cmin: mn,
          cmax: mx,
        },
        coordFor,
        byRootRep,
        selectedReps,
        repRanks
      ),
    ];
  }, [collectMetric, coordFor, byRootRep, selectedReps, repRanks]);

  const ecdFinalTraces = useMemo<Partial<ScatterData>[]>(() => {
    const vals = collectMetric((r) => r.ecd_ll_final);
    const [mn, mx] = extent(vals);
    return [
      buildConnectionsTrace(coordFor, byRootRep, selectedReps, repRanks),
      buildMarkersTrace(
        {
          title: "ecd_ll_final",
          values: (r) => r.ecd_ll_final,
          colorScale: REDS_SCALE, // light -> dark
          cmin: mn,
          cmax: mx,
        },
        coordFor,
        byRootRep,
        selectedReps,
        repRanks
      ),
    ];
  }, [collectMetric, coordFor, byRootRep, selectedReps, repRanks]);

  return (
    <div className="min-h-screen bg-black text-white p-6 max-w-[1400px] mx-auto space-y-4">
      <h1 className="text-2xl font-bold">Spiral Plots</h1>

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
                method === m ? "bg-white text-black" : "bg-gray-800 text-white hover:bg-gray-700"
              }`}
            >
              {m}
            </button>
          ))}
        </div>

        {/* Twist control */}
        <div className="flex items-center gap-2 ml-6">
          <label className="text-sm font-medium">Twist (°/rank)</label>
          <input
            type="range"
            min={0}
            max={60}
            step={1}
            value={twistDeg}
            onChange={(e) => setTwistDeg(Number(e.target.value))}
            className="w-40"
            aria-label="Twist degrees per rank"
          />
        </div>
        <input
          type="number"
          className="w-20 border rounded px-2 py-1 bg-black text-white border-gray-600"
          min={0}
          max={180}
          step={1}
          value={twistDeg}
          onChange={(e) => setTwistDeg(Number(e.target.value))}
        />

        {/* Resample */}
        <button
          type="button"
          onClick={onResample}
          className="ml-auto px-4 py-2 rounded bg-white text-black hover:bg-gray-100 disabled:opacity-50"
          disabled={!repsAll.length || loading}
        >
          Resample 5 reps
        </button>
      </div>

      {loading ? (
        <div className="p-4 border border-gray-700 rounded bg-black text-white">Loading…</div>
      ) : err ? (
        <div className="p-4 border border-gray-700 rounded text-red-400 bg-black">{err}</div>
      ) : !rows.length || !selectedReps.length ? (
        <div className="p-4 border border-gray-700 rounded bg-black text-white">No data</div>
      ) : (
        <>
          {/* Row 1: ll_init | ll_final */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <ChartCard title={`${method}: ll_init`} traces={initTraces} />
            <ChartCard title={`${method}: ll_final`} traces={finalTraces} />
          </div>

          {/* Row 2: ecd_ll_first | ecd_ll_final */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <ChartCard title={`${method}: ecd_ll_first`} traces={ecdFirstTraces} />
            <ChartCard title={`${method}: ecd_ll_final`} traces={ecdFinalTraces} />
          </div>
        </>
      )}
    </div>
  );

  // ---------- helpers ----------

  function extent(nums: number[]): [number, number] {
    if (!nums.length) return [0, 1];
    let mn = nums[0], mx = nums[0];
    for (const v of nums) {
      if (v < mn) mn = v;
      if (v > mx) mx = v;
    }
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
}

function buildConnectionsTrace(
  coordFor: (node: string, rank: number) => { x: number; y: number },
  byRootRep: Map<string, Map<number, SpiralRow>>,
  selectedReps: number[],
  repRanks: Map<number, number>
): Partial<ScatterData> {
  const xs: number[] = [];
  const ys: number[] = [];

  for (const node of NODES) {
    const perRep = byRootRep.get(node);
    if (!perRep) continue;

    const present = selectedReps
      .filter((rep) => perRep.has(rep))
      .map((rep) => ({ rep, rank: repRanks.get(rep) ?? 0 }))
      .sort((a, b) => a.rank - b.rank);

    if (present.length < 2) continue;

    for (const { rank } of present) {
      const { x, y } = coordFor(node, rank);
      xs.push(x);
      ys.push(y);
    }
    xs.push(NaN);
    ys.push(NaN);
  }

  return {
    type: "scatter",
    mode: "lines",
    x: xs,
    y: ys,
    line: { color: "rgba(255,255,255,0.35)", width: 1, dash: "dot" },
    hoverinfo: "skip",
    showlegend: false,
    name: "replicate-links",
  };
}

function buildMarkersTrace(
  cfg: {
    title: string;
    values: (r: SpiralRow) => number | null;
    colorScale: [number, string][];
    cmin: number;
    cmax: number;
  },
  coordFor: (node: string, rank: number) => { x: number; y: number },
  byRootRep: Map<string, Map<number, SpiralRow>>,
  selectedReps: number[],
  repRanks: Map<number, number>
): Partial<ScatterData> {
  const xs: number[] = [];
  const ys: number[] = [];
  const colors: number[] = [];
  const texts: string[] = [];

  for (const node of NODES) {
    const perRep = byRootRep.get(node);
    if (!perRep) continue;
    for (const rep of selectedReps) {
      const row = perRep.get(rep);
      if (!row) continue;
      const rank = repRanks.get(rep) ?? 0;
      const c = coordFor(node, rank);
      const val = cfg.values(row);
      if (val == null || !Number.isFinite(val)) continue;

      xs.push(c.x);
      ys.push(c.y);
      colors.push(val);
      texts.push(`${node}, rep ${rep}<br>${cfg.title}: ${val.toFixed(3)}`);
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
      colorscale: cfg.colorScale, // explicit light→dark
      reversescale: false,        // do NOT flip
      cmin: cfg.cmin,
      cmax: cfg.cmax,
      showscale: true,
      colorbar: {
        title: { text: cfg.title, font: { color: "white" } },
        tickfont: { color: "white" },
        thickness: 12,
        bgcolor: "black",
        outlinecolor: "white",
      },
    },
    name: cfg.title,
    showlegend: false,
  };
}

// Small card wrapper for each plot
function ChartCard({
  title,
  traces,
}: {
  title: string;
  traces: Partial<Data>[];
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
  };
  const config: Partial<Config> = { displayModeBar: false, responsive: true };

  return (
    <div className="bg-black border border-gray-700 rounded p-2">
      <Plot data={traces} layout={layout} config={config} style={{ width: "100%", height: "100%" }} />
    </div>
  );
}
