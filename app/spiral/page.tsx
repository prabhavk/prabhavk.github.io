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

// Color scales (sequential: light for small, dark for large)
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

/* -------------------- Choreography -------------------- */
type ChoreoBounce = {
  type: "bounce";
  name?: string;
  min: number;       // deg/rank
  max: number;       // deg/rank
  step: number;      // deg per tick
  everyMs: number;   // tick interval
  repeats: number;   // number of direction reversals to perform before advancing
  startAt?: "min" | "max"; // default "min"
};
type ChoreoHold = {
  type: "hold";
  name?: string;
  at: number;        // deg/rank
  durationMs: number;
};
type ChoreoSeg = ChoreoBounce | ChoreoHold;

/** Peacock-inspired choreographies */
const PEACOCK_CHOREOS: Record<
  "FanDisplay" | "CourtshipStrut" | "RainRipple" | "IridescentFlutter" | "RoyalGlide",
  ChoreoSeg[]
> = {
  FanDisplay: [
    { type: "hold",   name: "present",  at: 0,   durationMs: 500 },
    { type: "bounce", name: "fan-open", min: -4,  max: 4,  step: 1, everyMs: 260, repeats: 4, startAt: "min" },
    { type: "bounce", name: "display",  min: -12, max: 12, step: 2, everyMs: 200, repeats: 6, startAt: "min" },
    { type: "hold",   name: "proud",    at: 12,  durationMs: 700 },
    { type: "bounce", name: "sweep",    min: -16, max: 16, step: 3, everyMs: 160, repeats: 8, startAt: "min" },
    { type: "hold",   name: "still",    at: 0,   durationMs: 600 },
  ],
  CourtshipStrut: [
    { type: "bounce", name: "start",   min: -6,  max: 6,  step: 1, everyMs: 300, repeats: 3, startAt: "min" },
    { type: "bounce", name: "stride",  min: -10, max: 10, step: 2, everyMs: 220, repeats: 4, startAt: "min" },
    { type: "hold",   name: "beat",    at: -8,  durationMs: 500 },
    { type: "bounce", name: "parade",  min: -16, max: 16, step: 3, everyMs: 160, repeats: 6, startAt: "min" },
    { type: "hold",   name: "pose",    at: 0,   durationMs: 800 },
  ],
  RainRipple: [
    { type: "bounce", name: "drizzle", min: -3,  max: 3,  step: 1, everyMs: 180, repeats: 6, startAt: "min" },
    { type: "bounce", name: "splash",  min: -12, max: 12, step: 4, everyMs: 120, repeats: 4, startAt: "min" },
    { type: "hold",   name: "drop",    at: 0,   durationMs: 500 },
    { type: "bounce", name: "ripples", min: -8,  max: 8,  step: 2, everyMs: 160, repeats: 5, startAt: "min" },
    { type: "hold",   name: "shimmer", at: 6,   durationMs: 400 },
    { type: "bounce", name: "afterglow", min: -16, max: 16, step: 2, everyMs: 200, repeats: 4, startAt: "min" },
  ],
  IridescentFlutter: [
    { type: "bounce", name: "flicker-1", min: -2,  max: 2,  step: 1, everyMs: 90,  repeats: 8, startAt: "min" },
    { type: "hold",   name: "blink",     at: 0,   durationMs: 200 },
    { type: "bounce", name: "flicker-2", min: -5,  max: 5,  step: 1, everyMs: 110, repeats: 8, startAt: "min" },
    { type: "hold",   name: "blink",     at: 0,   durationMs: 200 },
    { type: "bounce", name: "gleam",     min: -10, max: 10, step: 3, everyMs: 130, repeats: 6, startAt: "min" },
    { type: "hold",   name: "flash",     at: 0,   durationMs: 300 },
    { type: "bounce", name: "flare",     min: -16, max: 16, step: 4, everyMs: 110, repeats: 4, startAt: "min" },
  ],
  RoyalGlide: [
    { type: "hold",   name: "poise",   at: -12, durationMs: 600 },
    { type: "bounce", name: "arc-L",   min: -14, max: 14, step: 2, everyMs: 340, repeats: 6, startAt: "min" },
    { type: "hold",   name: "crest",   at: 12,  durationMs: 600 },
    { type: "bounce", name: "arc-S",   min: -8,  max: 8,  step: 1, everyMs: 420, repeats: 6, startAt: "min" },
    { type: "hold",   name: "center",  at: 0,   durationMs: 800 },
    { type: "bounce", name: "royal",   min: -16, max: 16, step: 2, everyMs: 300, repeats: 8, startAt: "min" },
  ],
};

const CHOREO_BY_METHOD: Record<MethodName, ChoreoSeg[]> = {
  Dirichlet: PEACOCK_CHOREOS.CourtshipStrut,
  Parsimony: PEACOCK_CHOREOS.RainRipple,
  SSH:       PEACOCK_CHOREOS.FanDisplay,
};
/* ----------------------------------------------------- */

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

  // Derived choreography for current method
  const CHOREO = useMemo<ChoreoSeg[]>(() => CHOREO_BY_METHOD[method], [method]);

  // Twist & “dance” control (single toggle)
  const [twistDeg, setTwistDeg] = useState<number>(() => {
    const first = CHOREO_BY_METHOD.Dirichlet[0];
    if (first?.type === "bounce") return first.startAt === "max" ? first.max : first.min;
    if (first?.type === "hold") return first.at;
    return -16;
  });
  const [dancing, setDancing] = useState<boolean>(false);

  // Choreography runner refs
  const segIdxRef = useRef<number>(0);
  const timerRef = useRef<number | null>(null);
  const dirRef = useRef<1 | -1>(1);        // current direction for bounce
  const reversalsLeftRef = useRef<number>(0);

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

  // When switching methods: stop dancing & reset to starting pose
  useEffect(() => {
    setDancing(false);
    segIdxRef.current = 0;
    const first = CHOREO[0];
    if (first?.type === "bounce") setTwistDeg(first.startAt === "max" ? first.max : first.min);
    else if (first?.type === "hold") setTwistDeg(first.at);
  }, [CHOREO]);

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

    const dTwist = (Math.PI / 180) * twistDeg; // degrees -> radians per rank

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
      buildWebRingsTrace(coordFor, selectedReps.length),  
      buildMarkersTrace(
        { title: "ll_init", values: (r) => r.ll_init, colorScale: BLUES_SCALE, cmin: mn, cmax: mx },
        coordFor, byRootRep, selectedReps, repRanks
      ),
    ];
  }, [collectMetric, coordFor, byRootRep, selectedReps, repRanks]);

  const finalTraces = useMemo<Partial<ScatterData>[]>(() => {
    const vals = collectMetric((r) => r.ll_final);
    const [mn, mx] = extent(vals);
    return [
      buildConnectionsTrace(coordFor, byRootRep, selectedReps, repRanks),
      buildWebRingsTrace(coordFor, selectedReps.length),  
      buildMarkersTrace(
        { title: "ll_final", values: (r) => r.ll_final, colorScale: REDS_SCALE, cmin: mn, cmax: mx },
        coordFor, byRootRep, selectedReps, repRanks
      ),
    ];
  }, [collectMetric, coordFor, byRootRep, selectedReps, repRanks]);

  const ecdFirstTraces = useMemo<Partial<ScatterData>[]>(() => {
    const vals = collectMetric((r) => r.ecd_ll_first);
    const [mn, mx] = extent(vals);
    return [
      buildConnectionsTrace(coordFor, byRootRep, selectedReps, repRanks),
      buildWebRingsTrace(coordFor, selectedReps.length),  
      buildMarkersTrace(
        { title: "ecd_ll_first", values: (r) => r.ecd_ll_first, colorScale: BLUES_SCALE, cmin: mn, cmax: mx },
        coordFor, byRootRep, selectedReps, repRanks
      ),
    ];
  }, [collectMetric, coordFor, byRootRep, selectedReps, repRanks]);

  const ecdFinalTraces = useMemo<Partial<ScatterData>[]>(() => {
    const vals = collectMetric((r) => r.ecd_ll_final);
    const [mn, mx] = extent(vals);
    return [
      buildConnectionsTrace(coordFor, byRootRep, selectedReps, repRanks),
      buildWebRingsTrace(coordFor, selectedReps.length),  
      buildMarkersTrace(
        { title: "ecd_ll_final", values: (r) => r.ecd_ll_final, colorScale: REDS_SCALE, cmin: mn, cmax: mx },
        coordFor, byRootRep, selectedReps, repRanks
      ),
    ];
  }, [collectMetric, coordFor, byRootRep, selectedReps, repRanks]);

  /* -------------------- Choreography runner -------------------- */
  const stopTimer = useCallback(() => {
    if (timerRef.current) {
      window.clearInterval(timerRef.current);
      timerRef.current = null;
    }
  }, []);

  const startCurrentSegment = useCallback(() => {
    stopTimer();

    const seg = CHOREO[segIdxRef.current % CHOREO.length];
    if (!seg) return;

    if (seg.type === "hold") {
      setTwistDeg(seg.at);
      const t = window.setTimeout(() => {
        segIdxRef.current = (segIdxRef.current + 1) % CHOREO.length;
        startCurrentSegment();
      }, Math.max(50, seg.durationMs));
      // reuse timerRef for timeouts as well
      timerRef.current = t as unknown as number;
      return;
    }

    // bounce
    const startAt = seg.startAt === "max" ? seg.max : seg.min;
    setTwistDeg(startAt);
    dirRef.current = startAt === seg.max ? -1 : 1;
    reversalsLeftRef.current = Math.max(0, seg.repeats);

    const tick = () => {
      setTwistDeg((prev) => {
        let next = prev + dirRef.current * seg.step;
        if (next > seg.max) {
          next = seg.max;
          dirRef.current = -1;
          reversalsLeftRef.current -= 1;
        } else if (next < seg.min) {
          next = seg.min;
          dirRef.current = 1;
          reversalsLeftRef.current -= 1;
        }

        if (reversalsLeftRef.current <= 0) {
          stopTimer();
          segIdxRef.current = (segIdxRef.current + 1) % CHOREO.length;
          setTimeout(startCurrentSegment, 0);
        }
        return next;
      });
    };

    timerRef.current = window.setInterval(tick, Math.max(20, seg.everyMs));
  }, [CHOREO, stopTimer]);

  useEffect(() => {
    if (!dancing) {
      stopTimer();
      return;
    }
    const entry = CHOREO[segIdxRef.current % CHOREO.length];
    if (entry?.type === "bounce") setTwistDeg(entry.startAt === "max" ? entry.max : entry.min);
    else if (entry?.type === "hold") setTwistDeg(entry.at);
    startCurrentSegment();
    return () => stopTimer();
  }, [dancing, CHOREO, startCurrentSegment, stopTimer]);

  // Stop button handler: stop & reset to very first pose of current choreography
  const onStopDance = useCallback(() => {
    setDancing(false);
    segIdxRef.current = 0;
    const first = CHOREO[0];
    if (first?.type === "bounce") setTwistDeg(first.startAt === "max" ? first.max : first.min);
    else if (first?.type === "hold") setTwistDeg(first.at);
  }, [CHOREO]);

  /* ------------------------------------------------------------- */

  // Refs to each plotly graphDiv for export
  const initGDRef = useRef<PlotlyHTMLElement | null>(null);
  const finalGDRef = useRef<PlotlyHTMLElement | null>(null);
  const ecdFirstGDRef = useRef<PlotlyHTMLElement | null>(null);
  const ecdFinalGDRef = useRef<PlotlyHTMLElement | null>(null);

  // Export prefix (shared naming) — method appended in filename
  const exportPrefix = useMemo(
    () => buildExportPrefix({ job, twistDeg }),
    [job, twistDeg]
  );
  // const filenameFor = useCallback(
  //   (metric: string) => `${exportPrefix}__${method}__${metric}.png`,
  //   [exportPrefix, method]
  // );

  // Minimal toImage options type to avoid `any`
  type ToImageOpts = { format?: "png" | "svg" | "jpeg" | "webp"; height?: number; width?: number; scale?: number };

  const toImage = useCallback(async (gd: PlotlyHTMLElement, opts: ToImageOpts) => {
    const w = window as unknown as { Plotly?: { toImage: (el: PlotlyHTMLElement, o: ToImageOpts) => Promise<string> } };
    if (!w.Plotly?.toImage) throw new Error("Plotly.toImage not available");
    return w.Plotly.toImage(gd, opts);
  }, []);

  const downloadURI = (uri: string, name: string) => {
    const a = document.createElement("a");
    a.href = uri;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
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
    <div className="min-h-screen bg-black text-white p-6 max-w-[1400px] mx-auto space-y-4">
      <h1 className="text-2xl font-bold">Spiral Plots</h1>

      <div className="flex flex-wrap gap-3 items-end">
        <div className="text-sm text-gray-300">
          Job: <span className="font-mono">{job || "(none)"}</span>
        </div>

        {/* Method selector — switching disables dance */}
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

        {/* Status (read-only)
        <div className="ml-6 text-sm">
          Twist: <span className="font-mono">{twistDeg.toFixed(0)}°/rank</span>{" "}
          <span className="text-xs text-gray-500">
            segment {(segIdxRef.current % CHOREO.length) + 1}/{CHOREO.length}
          </span>
        </div> */}

        {/* Dance controls */}
        <div className="ml-auto flex gap-2">
          {dancing ? (
            <button
              type="button"
              onClick={onStopDance}
              className="px-4 py-2 rounded bg-red-500 text-white hover:bg-red-600"
              title="Stop and reset to starting pose"
            >
              Stop Dance
            </button>
          ) : (
            <button
              type="button"
              onClick={() => setDancing(true)}
              className="px-4 py-2 rounded bg-white text-black hover:bg-gray-100"
              title="Run the choreography"
            >
              Start Dance
            </button>
          )}

          {/* Resample (enabled when reps loaded & not loading) */}
          <button
            type="button"
            onClick={onResample}
            className="px-4 py-2 rounded bg-white text-black hover:bg-gray-100 disabled:opacity-50"
            disabled={!repsAll.length || loading}
          >
            Resample 5 reps
          </button>

          {/* Download all four plots for the current pose */}
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
          {/* Row 1: ll_init | ll_final */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <ChartCard title={`${method}: ll_init`} traces={initTraces} onReady={(gd) => (initGDRef.current = gd)} />
            <ChartCard title={`${method}: ll_final`} traces={finalTraces} onReady={(gd) => (finalGDRef.current = gd)} />
          </div>

          {/* Row 2: ecd_ll_first | ecd_ll_final */}
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

function buildWebRingsTrace(
  coordFor: (node: string, rank: number) => { x: number; y: number },
  maxRank: number
): Partial<ScatterData> {
  const xs: number[] = [];
  const ys: number[] = [];

  for (let rank = 0; rank < Math.max(0, maxRank); rank++) {
    if (NODES.length === 0) continue;

    // Walk around the ring
    for (let i = 0; i < NODES.length; i++) {
      const { x, y } = coordFor(NODES[i], rank);
      xs.push(x);
      ys.push(y);
    }

    // Close the loop by returning to the first node
    const first = coordFor(NODES[0], rank);
    xs.push(first.x);
    ys.push(first.y);

    // Break before next ring
    xs.push(NaN);
    ys.push(NaN);
  }

  return {
    type: "scatter",
    mode: "lines",
    x: xs,
    y: ys,
    line: { color: "rgba(255,255,255,0.20)", width: 1, dash: "dot" },
    hoverinfo: "skip",
    showlegend: false,
    name: "web-rings",
  };
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
    line: { color: "rgba(255,255,255,0.35)", width: 1, dash: "solid" },
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
      colorscale: cfg.colorScale,
      reversescale: false,
      cmin: cfg.cmin,
      cmax: cfg.cmax,
      showscale: true,
      colorbar: {
        title: { text: cfg.title, font: { color: "white" } },
        tickfont: { color: "white" },
        thickness: 12,
        tickformat: ".3f",
      },
    },
    name: cfg.title,
    showlegend: false,
  };
}
//
function ChartCard({
  title,
  traces,
  onReady,
}: {
  title: string;
  traces: Partial<Data>[];
  onReady?: (el: PlotlyHTMLElement | null) => void;
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
