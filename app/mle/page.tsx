// app/mle/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import dynamic from "next/dynamic";
import type { PlotParams } from "react-plotly.js";
// replace your plotly import with:
import type { Data, Layout, Config } from "plotly.js";


const Plot = dynamic<PlotParams>(() => import("react-plotly.js"), { ssr: false });

type EMStruct = {
  method?: string;
  rep?: number;
  ll_final?: number;
  root_name?: string;
  root_prob_init?: number[];
  root_prob_final?: number[];
  trans_prob_init?: unknown;   // unknown on purpose: we’ll normalize various shapes
  trans_prob_final?: unknown;  // "
};

type BestResp = {
  job_id: string;
  pars?: EMStruct | null;
  dirichlet?: EMStruct | null;
  ssh?: EMStruct | null;
};

type ApiErr = { error: string };

function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isApiErr(x: unknown): x is ApiErr {
  return isRecord(x) && typeof x.error === "string";
}
function isBestResp(x: unknown): x is BestResp {
  return isRecord(x) && typeof x.job_id === "string";
}

/* -------------------- math helpers (Beta pdf) -------------------- */
// Log-Gamma (Lanczos)
function logGamma(z: number): number {
  const p = [
    676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012,
    9.9843695780195716e-6, 1.5056327351493116e-7,
  ];
  if (z < 0.5) {
    return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - logGamma(1 - z);
  }
  z -= 1;
  let x = 0.99999999999980993;
  for (let i = 0; i < p.length; i++) x += p[i] / (z + i + 1);
  const t = z + p.length - 0.5;
  return 0.9189385332046727 + (z + 0.5) * Math.log(t) - t + Math.log(x); // 0.5*ln(2π)=0.918...
}
function logBeta(a: number, b: number): number {
  return logGamma(a) + logGamma(b) - logGamma(a + b);
}
function betaPdf(x: number, a: number, b: number): number {
  if (x <= 0 || x >= 1) return 0;
  const lp = (a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) - logBeta(a, b);
  return Math.exp(lp);
}
function linspace(n: number): number[] {
  return Array.from({ length: n }, (_, i) => i / (n - 1));
}
// Method-of-moments fit for Beta(a,b) from samples in (0,1)
function fitBetaParams(samples: number[]): { a: number; b: number } | null {
  const x = samples.filter((v) => Number.isFinite(v) && v > 0 && v < 1);
  if (x.length < 3) return null;
  const mean = x.reduce((s, v) => s + v, 0) / x.length;
  const var_ = x.reduce((s, v) => s + (v - mean) * (v - mean), 0) / (x.length - 1);

  // Guard against degenerate variance
  const eps = 1e-6;
  if (var_ <= eps || mean <= 0 || mean >= 1) return null;

  const common = (mean * (1 - mean)) / var_ - 1;
  const a = Math.max(eps, mean * common);
  const b = Math.max(eps, (1 - mean) * common);
  if (!Number.isFinite(a) || !Number.isFinite(b)) return null;
  return { a, b };
}

/* -------------------- transition-matrix normalization -------------------- */
/** Try to extract a 4x4 as rows from various possible JSON shapes */
function extractRows4x4(md: unknown): number[][] | null {
  // 1) Already an array of 4 arrays
  if (Array.isArray(md) && md.length === 4 && md.every((r) => Array.isArray(r) && r.length === 4)) {
    return md as number[][];
  }

  if (isRecord(md)) {
    // 2) { rows: number[][] } / { data: number[][] } / { m: number[][] }
    for (const key of ["rows", "data", "m"]) {
      const val = md[key];
      if (Array.isArray(val) && val.length === 4 && val.every((r) => Array.isArray(r) && r.length === 4)) {
        return val as number[][];
      }
    }

    // 3) { "0":[..], "1":[..], "2":[..], "3":[..] }
    const numKeys = ["0", "1", "2", "3"];
    if (numKeys.every((k) => Array.isArray(md[k]) && (md[k] as unknown[]).length === 4)) {
      return numKeys.map((k) => (md[k] as number[]).map(Number));
    }

    // 4) { r1:[..], r2:[..], r3:[..], r4:[..] }
    const rKeys = ["r1", "r2", "r3", "r4"];
    if (rKeys.every((k) => Array.isArray(md[k]) && (md[k] as unknown[]).length === 4)) {
      return rKeys.map((k) => (md[k] as number[]).map(Number));
    }

    // 5) { m11:.., m12:.., ... m44:.. }
    const hasM11 = Object.prototype.hasOwnProperty.call(md, "m11");
    const hasM44 = Object.prototype.hasOwnProperty.call(md, "m44");
    if (hasM11 && hasM44) {
      const rows: number[][] = [];
      for (let i = 1; i <= 4; i++) {
        const row: number[] = [];
        for (let j = 1; j <= 4; j++) {
          const key = `m${i}${j}`;
          const v = Number((md as Record<string, unknown>)[key]);
          row.push(Number.isFinite(v) ? v : NaN);
        }
        rows.push(row);
      }
      if (rows.every((r) => r.length === 4)) return rows;
    }
  }

  // 6) Flat array length 16 -> chunk
  if (Array.isArray(md) && md.length === 16) {
    const arr = (md as number[]).map(Number);
    const rows = [arr.slice(0, 4), arr.slice(4, 8), arr.slice(8, 12), arr.slice(12, 16)];
    return rows;
  }

  return null;
}

function collectDiagonalProbs(transFinal: unknown): number[] {
  const out: number[] = [];
  if (!isRecord(transFinal)) return out;

  for (const v of Object.values(transFinal)) {
    const rows = extractRows4x4(v);
    if (!rows) continue;
    // Diagonal p11, p22, p33, p44
    for (let i = 0; i < 4; i++) {
      const val = Number(rows[i][i]);
      if (Number.isFinite(val) && val >= 0 && val <= 1) out.push(val);
    }
  }
  return out;
}

/* -------------------- Page -------------------- */
const SELECTED_KEY = "emtr:selectedJobId";

export default function MLEinDPage() {
  const [job, setJob] = useState<string>("");
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);
  const [diagVals, setDiagVals] = useState<number[]>([]);

  // read selected job
  useEffect(() => {
    try {
      const saved = localStorage.getItem(SELECTED_KEY) ?? "";
      setJob(saved);
    } catch {
      /* ignore */
    }
  }, []);

  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Pick a job on the Precomputed Results page.");
      setDiagVals([]);
      return;
    }
    setLoading(true);
    setErr(null);
    try {
      const u = new URL("/api/best", window.location.origin);
      u.searchParams.set("job", job);
      const res = await fetch(u.toString(), { cache: "no-store" });
      const ct = res.headers.get("content-type") || "";
      if (!ct.includes("application/json")) {
        const body = await res.text().catch(() => "");
        throw new Error(`Expected JSON, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`);
      }
      const json: unknown = await res.json();
      if (!res.ok || isApiErr(json) || !isBestResp(json)) {
        const msg = isApiErr(json) ? json.error : `HTTP ${res.status}`;
        throw new Error(msg);
      }

      const diri = json.dirichlet;
      if (!diri || !diri.trans_prob_final) {
        setDiagVals([]);
        setErr("No EM_best_diri available for this job.");
      } else {
        const vals = collectDiagonalProbs(diri.trans_prob_final);
        setDiagVals(vals);
      }
    } catch (e: unknown) {
      setErr(e instanceof Error ? e.message : "Failed to load");
      setDiagVals([]);
    } finally {
      setLoading(false);
    }
  }, [job]);

  useEffect(() => {
    if (job) void load();
  }, [job, load]);

  // Fit Beta to the diagonal values for overlay
  const betaParams = useMemo(() => fitBetaParams(diagVals), [diagVals]);

  // Build traces: histogram (density) + Beta pdf line
  const traces = useMemo<Partial<Data>[]>(() => {
  const out: Partial<Data>[] = [];

  // histogram (with nbinsx) — cast through unknown to satisfy TS
  const hist = {
    type: "histogram" as const,
    x: diagVals,
    name: "Diagonal probabilities",
    histnorm: "probability density" as const,
    nbinsx: 30,
    marker: { color: "#1f77b4", line: { width: 1, color: "#fff" } },
    hovertemplate: "p<sub>ii</sub>=%{x:.3f}<br>density=%{y:.3f}<extra></extra>",
  };
  out.push(hist as unknown as Partial<Data>);

  // beta pdf overlay
  if (betaParams) {
    const xs = linspace(401);
    const ys = xs.map((x) => betaPdf(x, betaParams.a, betaParams.b));
    const line = {
      type: "scatter" as const,
      mode: "lines" as const,
      x: xs,
      y: ys,
      name: `Beta pdf (a=${betaParams.a.toFixed(2)}, b=${betaParams.b.toFixed(2)})`,
      line: { width: 3, color: "black" },
      hovertemplate: "x=%{x:.3f}<br>pdf=%{y:.3f}<extra></extra>",
    };
    out.push(line as unknown as Partial<Data>);
  }

  return out;
}, [diagVals, betaParams]);

  const layout: Partial<Layout> = {
    title: { text: "MLE in D · Diagonal transition probabilities (EM_best_diri)" },
    xaxis: { title: { text: "p(ii)" }, range: [0, 1] },
    yaxis: { title: { text: "Density" } },
    bargap: 0.05,
    margin: { l: 60, r: 20, t: 50, b: 50 },
    showlegend: true,
  };

  const config: Partial<Config> = { displayModeBar: false, responsive: true };

  return (
    <div className="p-6 max-w-[1200px] mx-auto">
      <h1 className="text-2xl font-bold mb-4">MLE in D</h1>

      <div className="flex flex-wrap gap-3 items-end mb-4">
        <div className="text-sm">
          Job: <span className="font-mono">{job || "(none)"}</span>
        </div>
        <button
          className="ml-auto px-4 py-2 rounded bg-black text-white disabled:opacity-50"
          onClick={() => void load()}
          disabled={!job || loading}
        >
          {loading ? "Loading…" : "Reload"}
        </button>
      </div>

      {loading ? (
        <div className="p-4 border rounded">Loading…</div>
      ) : err ? (
        <div className="p-4 border rounded text-red-600">{err}</div>
      ) : diagVals.length === 0 ? (
        <div className="p-4 border rounded">No diagonal probabilities found for EM_best_diri.</div>
      ) : (
        <div className="bg-white border rounded p-2">
          <Plot data={traces as Data[]} layout={layout} config={config} style={{ width: "100%", height: "100%" }} />
          <div className="text-sm text-gray-600 mt-2">
            Using histogram with <code>histnorm=&quot;probability density&quot;</code> so the Beta pdf overlays on the same scale.
          </div>
        </div>
      )}
    </div>
  );
}
