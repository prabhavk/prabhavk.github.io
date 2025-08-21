// app/spiral/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import dynamic from "next/dynamic";
import type { PlotParams } from "react-plotly.js";
import type { Data, Layout, Config, ScatterData } from "plotly.js";

const Plot = dynamic<PlotParams>(() => import("react-plotly.js"), { ssr: false });

type MethodName = "Parsimony" | "Dirichlet" | "SSH";

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

const BLUE_SCALE = "Blues";
const RED_SCALE = "Reds";

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
    (x["method"] === "Parsimony" || x["method"] === "Dirichlet" || x["method"] === "SSH") &&
    Array.isArray(x["rows"]) &&
    Array.isArray(x["reps"])
  );
}

export default function SpiralPage() {
  const [job, setJob] = useState<string>("");
  const [method, setMethod] = useState<MethodName>("Parsimony");
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
      const u = new URL("/api/flower", window.location.origin); // API path kept as /api/flower per your backend
      u.searchParams.set("job", job);
      u.searchParams.set("method", method);
      const res = await fetch(u.toString(), { cache: "no-store" });

      const ct = res.headers.get("content-type") || "";
      if (!ct.includes("application/json")) {
        const body = await res.text().catch(() => "");
        throw new Error(
          `Expected JSON from ${u.pathname}, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(
            0,
            160
          )}`
        );
      }

      const json: unknown = await res.json();
      if (!res.ok || isApiErr(json)) throw new Error(isApiErr(json) ? json.error : `HTTP ${res.status}`);
      if (!isApiSpiral(json)) throw new Error("Invalid /api/flower response");
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
    const r0 = 0.6; // base radius
    const dr = 0.35; // per-rank radial step
    const dTwist = Math.PI / 10; // twist per rank

    NODES.forEach((node, armIdx) => {
      const base = armIdx * dTheta;
      for (let k = 0; k < 5; k++) {
        const r = r0 + k * dr;
        const th = -(base + k * dTwist); // NEGATED => clockwise
        out.set(`${node}#${k}`, { x: r * Math.cos(th), y: r * Math.sin(th) });
      }
    });
    return out;
  }, []);

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
          colorScale: BLUE_SCALE,
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
          colorScale: RED_SCALE,
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
          colorScale: BLUE_SCALE,
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
          colorScale: RED_SCALE,
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
    <div className="p-6 max-w-[1400px] mx-auto space-y-4">
      <h1 className="text-2xl font-bold">Spiral Plots</h1>

      <div className="flex flex-wrap gap-3 items-end">
        <div className="text-sm text-gray-700">
          Job: <span className="font-mono">{job || "(none)"}</span>
        </div>

        {/* Method selector */}
        <div className="flex gap-2 ml-4">
          {(["Parsimony", "Dirichlet", "SSH"] as const).map((m) => (
            <button
              key={m}
              type="button"
              onClick={() => setMethod(m)}
              className={`px-3 py-2 rounded border ${
                method === m ? "bg-black text-white" : "bg-white hover:bg-gray-100"
              }`}
            >
              {m}
            </button>
          ))}
        </div>

        {/* Resample */}
        <button
          type="button"
          onClick={onResample}
          className="ml-auto px-4 py-2 rounded bg-black text-white disabled:opacity-50"
          disabled={!repsAll.length || loading}
        >
          Resample 5 reps
        </button>
      </div>

      {loading ? (
        <div className="p-4 border rounded">Loading…</div>
      ) : err ? (
        <div className="p-4 border rounded text-red-600">{err}</div>
      ) : !rows.length || !selectedReps.length ? (
        <div className="p-4 border rounded">No data</div>
      ) : (
        <>
          {/* Row 1: ll_init | ll_final */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <ChartCard title={`${method}: ll_init (before)`} traces={initTraces} />
            <ChartCard title={`${method}: ll_final (after)`} traces={finalTraces} />
          </div>

          {/* Row 2: ecd_ll_first | ecd_ll_final */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <ChartCard title={`${method}: ecd_ll_first (before)`} traces={ecdFirstTraces} />
            <ChartCard title={`${method}: ecd_ll_final (after)`} traces={ecdFinalTraces} />
          </div>
        </>
      )}
    </div>
  );

  // ---------- helpers ----------

  function extent(nums: number[]): [number, number] {
    if (!nums.length) return [0, 1];
    let mn = nums[0],
      mx = nums[0];
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

  // connect selected replicates for each node (in ascending rank)
  for (const node of NODES) {
    const perRep = byRootRep.get(node);
    if (!perRep) continue;

    const present = selectedReps
      .filter((rep) => perRep.has(rep))
      .map((rep) => ({ rep, rank: repRanks.get(rep) ?? 0 }))
      .sort((a, b) => a.rank - b.rank);

    if (present.length < 2) continue; // nothing to connect

    for (const { rank } of present) {
      const { x, y } = coordFor(node, rank);
      xs.push(x);
      ys.push(y);
    }
    // separator between nodes
    xs.push(NaN);
    ys.push(NaN);
  }

  return {
    type: "scatter",
    mode: "lines",
    x: xs,
    y: ys,
    line: { color: "rgba(0,0,0,0.35)", width: 1, dash: "dot" },
    hoverinfo: "skip",
    showlegend: false,
    name: "replicate-links",
  };
}

function buildMarkersTrace(
  cfg: {
    title: string;
    values: (r: SpiralRow) => number | null;
    colorScale: string;
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
      colorscale: cfg.colorScale,
      cmin: cfg.cmin,
      cmax: cfg.cmax,
      showscale: true,
      colorbar: { title: { text: cfg.title }, thickness: 12 },
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
    title: { text: title },
    xaxis: { visible: false, scaleanchor: "y" },
    yaxis: { visible: false },
    margin: { l: 20, r: 60, t: 40, b: 20 },
    height: 520,
    showlegend: false,
  };
  const config: Partial<Config> = { displayModeBar: false, responsive: true };

  return (
    <div className="bg-white border rounded p-2">
      <Plot data={traces} layout={layout} config={config} style={{ width: "100%", height: "100%" }} />
    </div>
  );
}
