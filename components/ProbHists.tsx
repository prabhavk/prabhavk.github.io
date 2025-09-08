// components/ProbHists.tsx
"use client";

import React, { useEffect, useMemo, useState } from "react";
import dynamic from "next/dynamic";
import type { PlotParams } from "react-plotly.js";
import type { PlotData, Layout, Config, Shape } from "plotly.js";

const Plot = dynamic<PlotParams>(() => import("react-plotly.js"), { ssr: false });

/* --------- colors --------- */
const COLOR_BLACK = "#000000";
const COLOR_GRAY_DENSITY = "#6b7280"; // neutral gray for density curve
const COLOR_CI = "#9CA3AF";           // slightly lighter gray for CI lines

/* --------- constants & layouts --------- */

const DNA = ["a", "c", "g", "t"] as const;

type HistogramXBins = { start: number | string; end: number | string; size: number | string };
const XBINS_0_1: HistogramXBins = { start: 0, end: 1, size: 0.02 }; // bin size

// Shared axis-frame styling (black frame around plotting area)
const axisFrame = {
  ticks: "outside" as const,
  showline: true,
  linecolor: "white",
  linewidth: 2,
  mirror: true as const,
  zeroline: false,
};

const smallLayout: Partial<Layout> = {
  margin: { l: 40, r: 10, t: 10, b: 35 },
  xaxis: { title: { text: "p" }, range: [0, 1], ...axisFrame },
  yaxis: { title: { text: "count" }, ...axisFrame },
  height: 240,
  showlegend: false,
  plot_bgcolor: "white",
  paper_bgcolor: "white",
};

const bigLayout: Partial<Layout> = {
  margin: { l: 40, r: 10, t: 10, b: 35 },
  xaxis: { title: { text: "p" }, range: [0, 1], ...axisFrame },
  yaxis: { title: { text: "count" }, ...axisFrame },
  height: 320,
  showlegend: false,
  plot_bgcolor: "white",
  paper_bgcolor: "white",
};

const baseConfig: Partial<Config> = { displayModeBar: false, responsive: true };

/* --------- labels --------- */

function label_trans_prob(row1: number, col1: number): string {
  const r = DNA[(row1 - 1) as 0 | 1 | 2 | 3] ?? "?";
  const c = DNA[(col1 - 1) as 0 | 1 | 2 | 3] ?? "?";
  return `p(${c}|${r})`;
}
function label_root_prob(idx1: number): string {
  const r = DNA[(idx1 - 1) as 0 | 1 | 2 | 3] ?? "?";
  return `p(${r})`;
}

/* --------- utils --------- */

function clamp01(v: number): number {
  if (!Number.isFinite(v)) return NaN;
  if (v < 0) return 0;
  if (v > 1) return 1;
  return v;
}
function isRec(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}

/** Accept many shapes and normalize to 4x4 matrix. */
function extractRows4x4(md: unknown): number[][] | null {
  if (Array.isArray(md) && md.length === 4 && md.every(r => Array.isArray(r) && r.length === 4)) {
    return (md as unknown[][]).map(r => r.map(n => clamp01(Number(n))));
  }
  if (isRec(md)) {
    for (const k of ["rows", "data", "m"] as const) {
      const v = md[k];
      if (Array.isArray(v) && v.length === 4 && (v as unknown[]).every(r => Array.isArray(r) && (r as unknown[]).length === 4)) {
        return (v as unknown[][]).map(r => r.map(n => clamp01(Number(n))));
      }
    }
    const keys01 = ["0", "1", "2", "3"] as const;
    if (keys01.every(k => Array.isArray(md[k]) && (md[k] as unknown[]).length === 4)) {
      return keys01.map(k => (md[k] as unknown[]).map(n => clamp01(Number(n))));
    }
    const keysR = ["r1", "r2", "r3", "r4"] as const;
    if (keysR.every(k => Array.isArray(md[k]) && (md[k] as unknown[]).length === 4)) {
      return keysR.map(k => (md[k] as unknown[]).map(n => clamp01(Number(n))));
    }
    const hasM11 = Object.prototype.hasOwnProperty.call(md, "m11");
    const hasM44 = Object.prototype.hasOwnProperty.call(md, "m44");
    if (hasM11 && hasM44) {
      const out: number[][] = [];
      for (let i = 1; i <= 4; i++) {
        const row: number[] = [];
        for (let j = 1; j <= 4; j++) {
          const v = clamp01(Number((md as Record<string, unknown>)[`m${i}${j}`]));
          row.push(v);
        }
        out.push(row);
      }
      return out;
    }
    const rec = md as Record<string, unknown>;
    const nr = Number((rec.nrow ?? rec.rows) as number);
    const nc = Number((rec.ncol ?? rec.cols) as number);
    const vals = (rec.values ?? rec.v ?? rec.flat) as unknown;
    if (nr === 4 && nc === 4 && Array.isArray(vals) && vals.length === 16) {
      const v = (vals as unknown[]).map(n => clamp01(Number(n)));
      return [v.slice(0,4), v.slice(4,8), v.slice(8,12), v.slice(12,16)];
    }
  }
  if (Array.isArray(md) && md.length === 16) {
    const v = (md as unknown[]).map(n => clamp01(Number(n)));
    return [v.slice(0,4), v.slice(4,8), v.slice(8,12), v.slice(12,16)];
  }
  return null;
}

/** Normalize map of matrices: Record<string, unknown> -> Record<string, number[][]> */
function normalizeTransMap(trans: Record<string, unknown> | null): Record<string, number[][]> {
  const out: Record<string, number[][]> = {};
  if (!trans || !isRec(trans)) return out;
  for (const [k, v] of Object.entries(trans)) {
    const m = extractRows4x4(v);
    if (m) out[k] = m;
  }
  return out;
}

/* --------- plotting helpers --------- */

function makeHistTrace(x: number[], name: string): Partial<PlotData> {
  return {
    type: "histogram",
    x,
    name,
    marker: { color: COLOR_BLACK, line: { width: 1, color: COLOR_BLACK } },
    opacity: 0.35, // black, slightly transparent for readability
    xbins: XBINS_0_1,
    hovertemplate: "p=%{x:.3f}<br>count=%{y}<extra></extra>",
  };
}

/* --------- Beta helpers (Dirichlet marginals) --------- */

function lnGamma(z: number): number {
  const g = 7;
  const C = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
  ];
  if (z < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - lnGamma(1 - z);
  z -= 1;
  let x = C[0];
  for (let i = 1; i < g + 2; i++) x += C[i] / (z + i);
  const t = z + g + 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}
function betaFunc(a: number, b: number): number {
  return Math.exp(lnGamma(a) + lnGamma(b) - lnGamma(a + b));
}
function betaPdf(x: number, a: number, b: number): number {
  if (x <= 0 || x >= 1) return 0;
  return Math.pow(x, a - 1) * Math.pow(1 - x, b - 1) / betaFunc(a, b);
}
function makeBetaLineTrace(a: number, b: number, totalCount: number, name: string): Partial<PlotData> {
  const binW = Number(XBINS_0_1.size) || 0.02;
  const N = Math.max(1, totalCount);
  const xs: number[] = [];
  const ys: number[] = [];
  const steps = 200;
  for (let i = 0; i <= steps; i++) {
    const x = i / steps;
    xs.push(x);
    const pdf = betaPdf(x, a, b);
    ys.push(pdf * binW * N);
  }
  return {
    type: "scatter",
    mode: "lines",
    x: xs,
    y: ys,
    name,
    line: { color: COLOR_GRAY_DENSITY, width: 2 },
    hovertemplate: "x=%{x:.3f}<br>exp.count/bin=%{y:.3f}<extra></extra>",
  };
}

/** Continued fraction for incomplete beta — Numerical Recipes style */
function betacf(a: number, b: number, x: number): number {
  const MAXIT = 200;
  const EPS = 3e-8;
  const FPMIN = 1e-30;

  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;
  let c = 1.0;
  let d = 1.0 - (qab * x) / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  let h = d;

  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;

    // even step
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    h *= d * c;

    // odd step
    aa = -((a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    const del = d * c;
    h *= del;

    if (Math.abs(del - 1.0) < EPS) break;
  }
  return h;
}

/** Regularized incomplete beta I_x(a,b) = P(X <= x) for Beta(a,b) */
function betaCdf(x: number, a: number, b: number): number {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(lnGamma(a + b) - lnGamma(a) - lnGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  if (x < (a + 1) / (a + b + 2)) {
    return (bt * betacf(a, b, x)) / a;
  } else {
    return 1 - (bt * betacf(b, a, 1 - x)) / b;
  }
}

/** Inverse CDF (quantile) via bisection */
function betaInvCDF(p: number, a: number, b: number): number {
  const pp = Math.min(1 - 1e-12, Math.max(1e-12, p));
  let lo = 0, hi = 1, mid = 0.5;
  for (let i = 0; i < 80; i++) {
    mid = 0.5 * (lo + hi);
    const c = betaCdf(mid, a, b);
    if (c > pp) hi = mid;
    else lo = mid;
    if (Math.abs(c - pp) < 1e-6) break;
  }
  return mid;
}

/** A vertical line that spans the whole plotting area (yref='paper') */
function makeVerticalShape(x: number, color = COLOR_BLACK, solid = true): Partial<Shape> {
  return {
    type: "line",
    x0: x,
    x1: x,
    y0: 0,
    y1: 1,
    xref: "x",
    yref: "paper",
    line: { color, width: 2, dash: solid ? "solid" : "dot" },
  };
}

/** Two dotted CI lines at given level for Beta(a,b) */
function makeBetaCIShapes(a: number, b: number, level = 0.95): Partial<Shape>[] {
  if (!(a > 0 && b > 0)) return [];
  const alpha = 1 - level;
  const lo = betaInvCDF(alpha / 2, a, b);
  const hi = betaInvCDF(1 - alpha / 2, a, b);
  return [makeVerticalShape(lo, COLOR_CI, false), makeVerticalShape(hi, COLOR_CI, false)];
}

/* --------- component --------- */

export function ProbHists({
  root,
  trans,
  rootKey,
  rootName,
  alphaPi,         // Dirichlet for root (length 4)
  alphaM,          // Dirichlet for transition rows (length 4): α1 (diag), α2=α3=α4 (off-diag)
  showBetaCI = false,
  betaCILevel = 0.95,
  rootHeaderRight,
  singleHeaderRight,
  aggHeaderRight,
}: {
  root: number[] | null;
  trans: Record<string, unknown> | null;
  rootKey?: string | null;
  rootName?: string | null;
  alphaPi?: number[] | null;
  alphaM?: number[] | null;
  showBetaCI?: boolean;
  betaCILevel?: number;
  rootHeaderRight?: React.ReactNode;
  singleHeaderRight?: React.ReactNode;
  aggHeaderRight?: React.ReactNode;
}) {
  // sanitize root
  const rootClean = useMemo(
    () => (Array.isArray(root) ? root.map(n => clamp01(Number(n))).filter(n => Number.isFinite(n)) : []),
    [root]
  );

  // normalize matrices & node keys
  const transMap = useMemo(() => normalizeTransMap(trans), [trans]);
  const nodeKeys = useMemo(() => {
    const keys = Object.keys(transMap).filter(k => (rootKey ? k !== rootKey : true));
    return keys.sort((a, b) => {
      const na = Number(a), nb = Number(b);
      const aNum = Number.isFinite(na), bNum = Number.isFinite(nb);
      if (aNum && bNum) return na - nb;
      if (aNum) return -1;
      if (bNum) return 1;
      return a.localeCompare(b);
    });
  }, [transMap, rootKey]);

  // selectors: row, col, node
  const [rowIdx, setRowIdx] = useState<number>(3);
  const [colIdx, setColIdx] = useState<number>(4);
  const [nodeKey, setNodeKey] = useState<string>("6");

  // default node
  useEffect(() => {
    if (!nodeKeys.length) setNodeKey("");
    else if (!nodeKey || !nodeKeys.includes(nodeKey)) setNodeKey(nodeKeys[0]);
  }, [nodeKeys, nodeKey]);

  // selected value (single node/entry)
  const selectedValue = useMemo(() => {
    if (!nodeKey || !(nodeKey in transMap)) return NaN;
    const r = Math.max(1, Math.min(4, rowIdx)) - 1;
    const c = Math.max(1, Math.min(4, colIdx)) - 1;
    const v = transMap[nodeKey][r]?.[c];
    return Number.isFinite(v) ? clamp01(Number(v)) : NaN;
  }, [transMap, nodeKey, rowIdx, colIdx]);

  // aggregated 16 arrays across nodes
  const transCells = useMemo(() => {
    const cells: number[][] = Array.from({ length: 16 }, () => []);
    for (const m of Object.values(transMap)) {
      for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
          const x = m[i][j];
          if (Number.isFinite(x)) cells[i * 4 + j].push(x);
        }
      }
    }
    return cells;
  }, [transMap]);

  const rootOK = rootClean.length === 4;
  const anyTrans = nodeKeys.length > 0;

  const hasAlphaPi = !!(alphaPi && alphaPi.length === 4 && alphaPi.every(a => Number.isFinite(a) && a > 0));
  const hasAlphaM  = !!(alphaM  && alphaM.length  === 4 && alphaM.every(a => Number.isFinite(a) && a > 0));
  const alpha0M    = hasAlphaM ? alphaM!.reduce((s, x) => s + x, 0) : 0;

  // shared wrapper class for black border around each histogram/plot
  const plotFrameCls = "border-2 border-black rounded";

  return (
    <div className="space-y-6">
      {/* Root panels */}
      <section>
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg text-gray-900 font-semibold">
            {rootName ? `Root distribution at ${rootName}` : "Root probability (final)"}
          </h2>
          <div>{rootHeaderRight ?? null}</div>
        </div>
        {rootOK ? (
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            {rootClean.map((val, i) => {
              const traces: Partial<PlotData>[] = [];
              const shapes: Partial<Shape>[] = [];

              if (hasAlphaPi) {
                const a_k = alphaPi![i];
                const a0  = alphaPi!.reduce((s, x) => s + x, 0);
                const b_k = a0 - a_k;
                traces.push(makeBetaLineTrace(a_k, b_k, 1, `Beta(${a_k.toFixed(1)}, ${b_k.toFixed(1)})`));
                if (showBetaCI) shapes.push(...makeBetaCIShapes(a_k, b_k, betaCILevel));
              }
              if (Number.isFinite(val)) shapes.push(makeVerticalShape(val));

              return (
                <div key={`root-${i}`} className="rounded p-5 border-2 border-black">
                  <div className={plotFrameCls}>
                    <Plot
                      data={traces}
                      layout={{
                        ...smallLayout,
                        bargap: 0.05,
                        shapes,
                        annotations: [
                          {
                            text: label_root_prob(i + 1) + (rootName ? ` @ root ${rootName}` : ""),
                            xref: "paper",
                            yref: "paper",
                            x: 0.5, y: 0.95,
                            showarrow: false,
                            font: { size: 14 },
                            align: "center",
                          },
                        ],
                      }}
                      config={baseConfig}
                      style={{ width: "100%", height: "100%" }}
                    />
                  </div>
                  <div className="text-xs text-gray-900 mt-0">
                    probability: <span className="font-mono">{val.toFixed(4)}</span>
                  </div>
                </div>
              );
            })}
          </div>
        ) : (
          <div className="text-sm text-gray-600">No root_prob_final (expected 4 values).</div>
        )}
      </section>

      {/* Single node-specific entry */}
      <section>
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg text-black font-semibold">Transition matrix entry (per node)</h2>
          <div>{singleHeaderRight ?? null}</div>
        </div>
        {!anyTrans ? (
          <div className="text-sm text-gray-600">No transition matrices available.</div>
        ) : (
          <div className="border-2 border-black rounded p-3 space-y-3">
            <div className="flex flex-wrap items-center gap-3">
              <label className="text-black text-sm">
                Child Node&nbsp;
                <select
                  className="border rounded px-2 py-1 text-black text-sm"
                  value={nodeKey}
                  onChange={(e) => setNodeKey(e.target.value)}
                >
                  {nodeKeys.map(k => (
                    <option key={k} value={k}>{k}</option>
                  ))}
                </select>
              </label>
              <label className="text-sm">
                Row&nbsp;
                <select
                  className="border rounded px-2 py-1 text-black text-sm"
                  value={rowIdx}
                  onChange={(e) => setRowIdx(Number(e.target.value))}
                >
                  {[1, 2, 3, 4].map(n => <option key={n} value={n}>{n}</option>)}
                </select>
              </label>
              <label className="text-sm">
                Col&nbsp;
                <select
                  className="border rounded px-2 py-1 text-black text-sm"
                  value={colIdx}
                  onChange={(e) => setColIdx(Number(e.target.value))}
                >
                  {[1, 2, 3, 4].map(n => <option key={n} value={n}>{n}</option>)}
                </select>
              </label>
              {rootKey ? (
                <div className="text-xs text-black">Root excluded: <code>{rootKey}</code></div>
              ) : null}
            </div>

            {(() => {
              const traces: Partial<PlotData>[] = [];
              const shapes: Partial<Shape>[] = [];
              if (hasAlphaM) {
                const isDiag = rowIdx === colIdx;
                const a = isDiag ? alphaM![0] : alphaM![1];
                const b = alpha0M - a;
                traces.push(makeBetaLineTrace(a, b, 1, `Beta(${a.toFixed(1)}, ${b.toFixed(1)})`));
                if (showBetaCI) shapes.push(...makeBetaCIShapes(a, b, betaCILevel));
              }
              if (Number.isFinite(selectedValue)) shapes.push(makeVerticalShape(selectedValue));

              return (
                <div className={plotFrameCls}>
                  <Plot
                    data={traces}
                    layout={{
                      ...bigLayout,
                      bargap: 0.05,
                      shapes,
                      annotations: [
                        {
                          text: `${label_trans_prob(rowIdx, colIdx)} for child node ${nodeKey || "?"}`,
                          xref: "paper",
                          yref: "paper",
                          x: 0.5, y: 0.95,
                          showarrow: false,
                          font: { size: 14 },
                          align: "center",
                        },
                      ],
                    }}
                    config={baseConfig}
                    style={{ width: "100%", height: "100%" }}
                  />
                </div>
              );
            })()}

            <div className="text-xs text-gray-900">
              {Number.isFinite(selectedValue)
                ? <>Selected probability:&nbsp;<span className="font-mono">{selectedValue.toFixed(4)}</span></>
                : <>No value available for this selection.</>}
            </div>
          </div>
        )}
      </section>

      {/* Aggregated across matrices */}
      <section>
        <div className="flex items-center justify-between mb-2">
          <h2 className="text-lg text-black font-semibold">Transition matrix entries (aggregated)</h2>
          <div>{aggHeaderRight ?? null}</div>
        </div>
        {transCells.some(arr => arr.length) ? (
          <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-4 gap-4">
            {transCells.map((arr, idx) => {
              const i = Math.floor(idx / 4) + 1;
              const j = (idx % 4) + 1;

              const traces: Partial<PlotData>[] = [ makeHistTrace(arr, `${label_trans_prob(i, j)}`) ];
              const shapes: Partial<Shape>[] = [];

              if (hasAlphaM) {
                const isDiag = i === j;
                const a = isDiag ? alphaM![0] : alphaM![1];
                const b = alpha0M - a;
                traces.push(makeBetaLineTrace(a, b, arr.length, `Beta(${a.toFixed(1)}, ${b.toFixed(1)})`));
                if (showBetaCI) shapes.push(...makeBetaCIShapes(a, b, betaCILevel));
              }

              return (
                <div key={`m-${idx}`} className="rounded p-2 border-2 border-black">
                  <div className={plotFrameCls}>
                    <Plot
                      data={traces}
                      layout={{
                        ...smallLayout,
                        bargap: 0.05,
                        shapes,
                        annotations: [
                          {
                            text: `${label_trans_prob(i, j)} (n=${arr.length})`,
                            xref: "paper",
                            yref: "paper",
                            x: 0.5, y: 0.95,
                            showarrow: false,
                            font: { size: 13 },
                            align: "center",
                          },
                        ],
                      }}
                      config={baseConfig}
                      style={{ width: "100%", height: "100%" }}
                    />
                  </div>
                </div>
              );
            })}
          </div>
        ) : (
          <div className="text-sm text-gray-900">
            No transition probabilities found in <code>trans_prob_final</code>.
          </div>
        )}
      </section>
    </div>
  );
}
