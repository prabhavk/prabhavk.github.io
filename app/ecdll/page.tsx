// app/ecdll/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import dynamic from "next/dynamic";
import type { Layout, PlotData } from "plotly.js";
import type { PlotParams } from "react-plotly.js";

const Plot = dynamic<PlotParams>(() => import("react-plotly.js"), { ssr: false });

/* ---------------- Types ---------------- */

type MethodKey = "dirichlet" | "parsimony" | "ssh";

const METHOD_LABEL: Record<MethodKey, string> = {
  dirichlet: "Dirichlet",
  parsimony: "Parsimony",
  ssh: "SSH",
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
  rep: number | null;
  points: Pt[];
  summary: SummaryRow | null;
};

type ApiErr = { ok: false; error: string };

/* ---------------- Helpers ---------------- */

/** Parse server payload into Pt[]. Accepts [[iter, ll], ...] or [{iter,ll}, ...]. */
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

/**
 * Gap^alpha transform for ECD-LL curves (composite plot only).
 * - Let best = max(ll). gap = best - ll >= 0.
 * - y' = (gap / maxGap)^alpha in [0,1]. (0 best, 1 worst/earliest)
 */
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

/**
 * Box–Cox transform (shared positive shift across all series)
 */
function boxCoxTransformThree(
  series: { name: string; pts: Pt[] | null }[],
  lambda: number,
  eps = 1e-6
): Record<string, Pt[]> {
  const allY: number[] = [];
  for (const s of series) {
    if (!s.pts) continue;
    for (const [, y] of s.pts) if (Number.isFinite(y)) allY.push(y);
  }
  if (allY.length === 0) return Object.fromEntries(series.map((s) => [s.name, []]));
  const minY = Math.min(...allY);
  const shift = minY <= 0 ? -minY + eps : 0;

  const out: Record<string, Pt[]> = {};
  for (const s of series) {
    if (!s.pts || !s.pts.length) {
      out[s.name] = [];
      continue;
    }
    const arr: Pt[] = s.pts.map(([it, y]) => {
      const yp = y + shift; // positive
      const yt = lambda === 0 ? Math.log(yp) : (Math.pow(yp, lambda) - 1) / lambda;
      return [it, yt];
    });
    out[s.name] = arr;
  }
  return out;
}

/** Build a monochrome (black) line trace with chosen dash style */
function monoTrace(name: string, pts: Pt[] | null, dash: "solid" | "dash" | "dot"): Partial<PlotData> | null {
  if (!pts || !pts.length) return null;
  return {
    type: "scatter",
    mode: "lines",
    name,
    x: pts.map((p) => p[0]),
    y: pts.map((p) => p[1]),
    line: { color: "#000000", dash, width: 2 },
  };
}

/* ---------------- Page ---------------- */

export default function EcdllPage() {
  const [job, setJob] = useState<string>("");
  const [err, setErr] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  // available reps and selected rep
  const [reps, setReps] = useState<number[]>([]);
  const [rep, setRep] = useState<number | null>(null);

  // Raw points per method (selected rep)
  const [dirPts, setDirPts] = useState<Pt[] | null>(null);
  const [parPts, setParPts] = useState<Pt[] | null>(null);
  const [sshPts, setSshPts] = useState<Pt[] | null>(null);

  // Summaries per method (selected rep)
  const [_dirSum, setDirSum] = useState<SummaryRow | null>(null);
  const [_parSum, setParSum] = useState<SummaryRow | null>(null);
  const [_sshSum, setSshSum] = useState<SummaryRow | null>(null);

  // Tuning for composite transform & Box–Cox
  const [alpha, setAlpha] = useState<number>(0.5);
  const [lambda, setLambda] = useState<number>(0.5); // Box–Cox λ

  // Refs to Plotly graph divs for downloads
  const graphRefs = useRef<Record<string, HTMLDivElement | null>>({
    boxcox: null,
    composite: null,
    dir: null,
    par: null,
    ssh: null,
  });

  // Read selected job id
  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId");
      if (saved) setJob(saved);
    } catch {
      /* ignore */
    }
  }, []);

  // Fetch list of reps for this job
  const loadReps = useCallback(
    async (jobId: string) => {
      const u = new URL("/api/ecdll/reps", window.location.origin);
      u.searchParams.set("job_id", jobId);
      const res = await fetch(u.toString(), { cache: "no-store" });
      const j = (await res.json()) as { ok: boolean; reps?: number[]; error?: string };
      if (!res.ok || !j.ok) throw new Error(j.error || `HTTP ${res.status}`);
      const list = j.reps ?? [];
      setReps(list);
      if (list.length && (rep === null || !list.includes(rep))) {
        setRep(list[0]); // default to first available rep
      }
    },
    [rep]
  );

  useEffect(() => {
    if (job) {
      loadReps(job).catch((e) => setErr(e instanceof Error ? e.message : "Failed to load repetitions"));
    }
  }, [job, loadReps]);

  const fetchOne = useCallback(
    async (method: MethodKey, repSel: number | null): Promise<{ pts: Pt[]; sum: SummaryRow | null }> => {
      const u = new URL("/api/ecdll", window.location.origin);
      u.searchParams.set("job_id", job);
      u.searchParams.set("method", method);
      if (repSel !== null) u.searchParams.set("rep", String(repSel));

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
        const msg = "error" in j ? j.error : `HTTP ${res.status}`;
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
      setDirPts(null); setParPts(null); setSshPts(null);
      setDirSum(null); setParSum(null); setSshSum(null);
      return;
    }
    if (rep === null) {
      setErr("No repetition selected.");
      return;
    }
    setErr(null);
    setLoading(true);
    try {
      const [d, p, s] = await Promise.allSettled([
        fetchOne("dirichlet", rep),
        fetchOne("parsimony", rep),
        fetchOne("ssh", rep),
      ]);

      setDirPts(d.status === "fulfilled" ? d.value.pts : []);
      setParPts(p.status === "fulfilled" ? p.value.pts : []);
      setSshPts(s.status === "fulfilled" ? s.value.pts : []);

      setDirSum(d.status === "fulfilled" ? d.value.sum : null);
      setParSum(p.status === "fulfilled" ? p.value.sum : null);
      setSshSum(s.status === "fulfilled" ? s.value.sum : null);

      const errs: string[] = [];
      if (d.status === "rejected") errs.push(d.reason?.message || "Dirichlet fetch failed");
      if (p.status === "rejected") errs.push(p.reason?.message || "Parsimony fetch failed");
      if (s.status === "rejected") errs.push(s.reason?.message || "SSH fetch failed");
      if (errs.length) setErr(errs.join(" | "));
    } catch (e) {
      setErr(e instanceof Error ? e.message : "Failed to load ECD-LL data");
      setDirPts([]); setParPts([]); setSshPts([]);
      setDirSum(null); setParSum(null); setSshSum(null);
    } finally {
      setLoading(false);
    }
  }, [job, rep, fetchOne]);

  // re-load when rep changes
  useEffect(() => {
    if (job && rep !== null) void load();
  }, [job, rep, load]);

  /* -------- Layouts (integer x ticks everywhere) -------- */

  const baseLayout: Partial<Layout> = useMemo(
    () => ({
      margin: { l: 60, r: 20, t: 30, b: 50 },
      xaxis: {
        title: { text: "iteration" },
        tickmode: "linear",
        dtick: 1,
        tickformat: "d",
      },
      yaxis: { title: { text: "ECD-LL" } },
      height: 320,
      showlegend: true,
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

  const boxCoxLayout: Partial<Layout> = useMemo(
    () => ({
      ...baseLayout,
      yaxis: { title: { text: `Box–Cox (λ=${lambda.toFixed(2)})` } },
      height: 360,
    }),
    [baseLayout, lambda]
  );

  /* -------- Build traces (monochrome styles) -------- */

  const compositeTraces = useMemo(() => {
    const traces: Partial<PlotData>[] = [];
    const push = (name: string, pts: Pt[] | null, dash: "solid" | "dash" | "dot") => {
      if (!pts || !pts.length) return;
      const t = transformGapPower(pts, alpha);
      traces.push({
        type: "scatter",
        mode: "lines",
        name,
        x: t.map((p) => p[0]),
        y: t.map((p) => p[1]),
        line: { color: "#000000", dash, width: 2 },
      });
    };
    push("Dirichlet", dirPts, "solid");
    push("Parsimony", parPts, "dash");
    push("SSH", sshPts, "dot");
    return traces;
  }, [dirPts, parPts, sshPts, alpha]);

  const boxCoxTraces = useMemo(() => {
    const transformed = boxCoxTransformThree(
      [
        { name: "Dirichlet", pts: dirPts },
        { name: "Parsimony", pts: parPts },
        { name: "SSH", pts: sshPts },
      ],
      lambda
    );
    const ds: Partial<PlotData>[] = [];
    const add = (name: string, dash: "solid" | "dash" | "dot") => {
      const pts = transformed[name] || [];
      if (!pts.length) return;
      ds.push({
        type: "scatter",
        mode: "lines",
        name,
        x: pts.map((p) => p[0]),
        y: pts.map((p) => p[1]),
        line: { color: "#000000", dash, width: 2 },
      });
    };
    add("Dirichlet", "solid");
    add("Parsimony", "dash");
    add("SSH", "dot");
    return ds;
  }, [dirPts, parPts, sshPts, lambda]);

  const dirTrace = useMemo(() => monoTrace("Dirichlet", dirPts, "solid"), [dirPts]);
  const parTrace = useMemo(() => monoTrace("Parsimony", parPts, "dash"), [parPts]);
  const sshTrace = useMemo(() => monoTrace("SSH", sshPts, "dot"), [sshPts]);

  /* -------- Download helpers (no 'plotly.js-dist-min' dependency) -------- */

  type PlotKey = "boxcox" | "composite" | "dir" | "par" | "ssh";

interface PlotlyLike {
  downloadImage: (
    gd: HTMLElement,
    opts: {
      format?: "png" | "jpeg" | "svg" | "webp";
      width?: number;
      height?: number;
      filename?: string;
    }
  ) => Promise<string>;
}

function isRecord(v: unknown): v is Record<string, unknown> {
  return typeof v === "object" && v !== null;
}

function isPlotlyLike(v: unknown): v is PlotlyLike {
  return isRecord(v) && typeof v["downloadImage"] === "function";
}

async function resolvePlotlyFrom(gd: HTMLDivElement): Promise<PlotlyLike> {
  // 1) Plotly already on window (UMD path)
  const win = (gd.ownerDocument?.defaultView ?? null) as (Window & { Plotly?: unknown }) | null;
  if (win?.Plotly && isPlotlyLike(win.Plotly)) return win.Plotly;

  // 2) Bundled dist (avoids Node polyfills)
  const mod = (await import("plotly.js-dist")) as unknown;
  const candidate =
    isRecord(mod) && "default" in mod ? ((mod as { default: unknown }).default as unknown) : mod;
  if (isPlotlyLike(candidate)) return candidate;

  throw new Error("Plotly not found. Install 'plotly.js-dist' or expose window.Plotly.");
}

const downloadFigure = React.useCallback(
  async (key: PlotKey, filename: string, w = 1600, h = 900) => {
    const gd = graphRefs.current[key];
    if (!gd) return;
    const Plotly = await resolvePlotlyFrom(gd);
    await Plotly.downloadImage(gd, { format: "png", width: w, height: h, filename });
  },
  []
);

const downloadAll = React.useCallback(async () => {
  const tasks: Array<Promise<void>> = [];
  const base = `ecdll_job-${job || "unknown"}_rep-${rep ?? "na"}`;

  if (graphRefs.current.boxcox)
    tasks.push(downloadFigure("boxcox", `${base}_boxcox_lambda-${lambda.toFixed(2)}`));
  if (graphRefs.current.composite)
    tasks.push(downloadFigure("composite", `${base}_composite_alpha-${alpha.toFixed(2)}`));
  if (graphRefs.current.dir) tasks.push(downloadFigure("dir", `${base}_dirichlet`));
  if (graphRefs.current.par) tasks.push(downloadFigure("par", `${base}_parsimony`));
  if (graphRefs.current.ssh) tasks.push(downloadFigure("ssh", `${base}_ssh`));

  // run sequentially to avoid overwhelming the browser
  for (const t of tasks) await t;
  }, [alpha, lambda, job, rep, downloadFigure]);

  /* ---------------- Render ---------------- */

  return (
    <div className="p-6 max-w-[1400px] mx-auto space-y-5">
      <div className="flex items-end gap-3">
        <h1 className="text-2xl font-bold">Expected Complete-Data LL Convergence</h1>
        <div className="text-sm ml-2">
          Job: <span className="font-mono">{job || "(none)"} </span>
        </div>

        {/* Repetition selector */}
        <label className="ml-6 text-sm flex items-center gap-2">
          Repetition:
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
          disabled={!job || rep === null || loading}
        >
          {loading ? "Loading…" : "Reload"}
        </button>

        {/* Download ALL */}
        <button
          type="button"
          onClick={() => void downloadAll()}
          className="px-3 py-1.5 rounded bg-black text-white disabled:opacity-50 text-sm"
          disabled={loading}
          title="Download all figures as PNG"
        >
          Download all PNGs
        </button>
      </div>

      {/* -------- Box–Cox overlay plot -------- */}
      <section className="border rounded p-3">
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg font-semibold">Box–Cox transform (overlay)</h2>
          <div className="flex items-center gap-3 text-sm">
            <label className="flex items-center gap-2">
              λ:
              <input
                type="range"
                min={0}
                max={1}
                step={0.05}
                value={lambda}
                onChange={(e) => setLambda(Number(e.target.value))}
              />
              <span className="text-gray-600">{lambda.toFixed(2)}</span>
            </label>
            <button
              type="button"
              onClick={() =>
                downloadFigure(
                  "boxcox",
                  `ecdll_job-${job || "unknown"}_rep-${rep ?? "na"}_boxcox_lambda-${lambda.toFixed(2)}`
                )
              }
              className="px-2 py-1 rounded bg-black text-white"
            >
              Download PNG
            </button>
          </div>
        </div>
        {boxCoxTraces.length ? (
          <Plot
            data={boxCoxTraces}
            layout={boxCoxLayout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
            onInitialized={(_, gd) => { graphRefs.current.boxcox = gd as HTMLDivElement; }}
            onUpdate={(_, gd) => { graphRefs.current.boxcox = gd as HTMLDivElement; }}
          />
        ) : (
          <div className="text-sm text-gray-600">No data available for Box–Cox plot.</div>
        )}
      </section>

      {/* Alpha control for composite transform */}
      <div className="flex items-center gap-3">
        <label className="text-sm">Composite transform α (0.2–0.9)</label>
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
      <section className="border rounded p-3">
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg font-semibold">Composite (Gap^α transform)</h2>
          <button
            type="button"
            onClick={() =>
              downloadFigure(
                "composite",
                `ecdll_job-${job || "unknown"}_rep-${rep ?? "na"}_composite_alpha-${alpha.toFixed(2)}`
              )
            }
            className="px-2 py-1 rounded bg-black text-white"
          >
            Download PNG
          </button>
        </div>
        {compositeTraces.length ? (
          <Plot
            data={compositeTraces}
            layout={compositeLayout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
            onInitialized={(_, gd) => { graphRefs.current.composite = gd as HTMLDivElement; }}
            onUpdate={(_, gd) => { graphRefs.current.composite = gd as HTMLDivElement; }}
          />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data found for any method.</div>
        )}
      </section>

      {/* Individual plots */}
      <section className="border rounded p-3">
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg font-semibold">Dirichlet</h2>
          <button
            type="button"
            onClick={() =>
              downloadFigure("dir", `ecdll_job-${job || "unknown"}_rep-${rep ?? "na"}_dirichlet`)
            }
            className="px-2 py-1 rounded bg-black text-white"
          >
            Download PNG
          </button>
        </div>
        {dirTrace ? (
          <Plot
            data={[dirTrace as Partial<PlotData>]}
            layout={baseLayout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
            onInitialized={(_, gd) => { graphRefs.current.dir = gd as HTMLDivElement; }}
            onUpdate={(_, gd) => { graphRefs.current.dir = gd as HTMLDivElement; }}
          />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      <section className="border rounded p-3">
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg font-semibold">Parsimony</h2>
          <button
            type="button"
            onClick={() =>
              downloadFigure("par", `ecdll_job-${job || "unknown"}_rep-${rep ?? "na"}_parsimony`)
            }
            className="px-2 py-1 rounded bg-black text-white"
          >
            Download PNG
          </button>
        </div>
        {parTrace ? (
          <Plot
            data={[parTrace as Partial<PlotData>]}
            layout={baseLayout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
            onInitialized={(_, gd) => { graphRefs.current.par = gd as HTMLDivElement; }}
            onUpdate={(_, gd) => { graphRefs.current.par = gd as HTMLDivElement; }}
          />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      <section className="border rounded p-3">
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg font-semibold">SSH</h2>
          <button
            type="button"
            onClick={() =>
              downloadFigure("ssh", `ecdll_job-${job || "unknown"}_rep-${rep ?? "na"}_ssh`)
            }
            className="px-2 py-1 rounded bg-black text-white"
          >
            Download PNG
          </button>
        </div>
        {sshTrace ? (
          <Plot
            data={[sshTrace as Partial<PlotData>]}
            layout={baseLayout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
            onInitialized={(_, gd) => { graphRefs.current.ssh = gd as HTMLDivElement; }}
            onUpdate={(_, gd) => { graphRefs.current.ssh = gd as HTMLDivElement; }}
          />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      {err ? <div className="p-3 border rounded text-red-600">{err}</div> : null}
    </div>
  );
}
