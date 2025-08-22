// components/ProbHists.tsx
"use client";

import React, { useMemo, useState } from "react";
import dynamic from "next/dynamic";
import type { PlotParams } from "react-plotly.js";
import type { PlotData, Layout, Config } from "plotly.js";

const Plot = dynamic<PlotParams>(() => import("react-plotly.js"), { ssr: false });

/* --------- layouts & config --------- */

const smallLayout: Partial<Layout> = {
  margin: { l: 40, r: 10, t: 30, b: 35 },
  xaxis: { title: { text: "p" }, range: [0, 1] },
  yaxis: { title: { text: "count" } },
  height: 240,
};

const bigLayout: Partial<Layout> = {
  margin: { l: 40, r: 10, t: 30, b: 35 },
  xaxis: { title: { text: "p" }, range: [0, 1] },
  yaxis: { title: { text: "count" } },
  height: 320,
};

const baseConfig: Partial<Config> = { displayModeBar: false, responsive: true };

type HistogramXBins = { start: number | string; end: number | string; size: number | string };
const XBINS_0_1: HistogramXBins = { start: 0, end: 1, size: 0.02 };

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

/** Normalize a variety of shapes into a 4x4 (row-major) matrix. */
function extractRows4x4(md: unknown): number[][] | null {
  // [[...],[...],[...],[...]]
  if (Array.isArray(md) && md.length === 4 && md.every(r => Array.isArray(r) && r.length === 4)) {
    return (md as unknown[][]).map(r => r.map(n => clamp01(Number(n))));
  }
  if (isRec(md)) {
    // { rows|data|m: number[][] }
    for (const k of ["rows", "data", "m"] as const) {
      const v = md[k];
      if (Array.isArray(v) && v.length === 4 && (v as unknown[]).every(r => Array.isArray(r) && (r as unknown[]).length === 4)) {
        return (v as unknown[][]).map(r => r.map(n => clamp01(Number(n))));
      }
    }
    // { "0":[..], "1":[..], "2":[..], "3":[..] }
    const keys01 = ["0", "1", "2", "3"] as const;
    if (keys01.every(k => Array.isArray(md[k]) && (md[k] as unknown[]).length === 4)) {
      return keys01.map(k => (md[k] as unknown[]).map(n => clamp01(Number(n))));
    }
    // { r1:[..], r2:[..], r3:[..], r4:[..] }
    const keysR = ["r1", "r2", "r3", "r4"] as const;
    if (keysR.every(k => Array.isArray(md[k]) && (md[k] as unknown[]).length === 4)) {
      return keysR.map(k => (md[k] as unknown[]).map(n => clamp01(Number(n))));
    }
    // { m11..m44 }
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
    // { nrow:4, ncol:4, values:[16] } (row-major)
    const rec = md as Record<string, unknown>;
    const nr = Number((rec.nrow ?? rec.rows) as number);
    const nc = Number((rec.ncol ?? rec.cols) as number);
    const vals = (rec.values ?? rec.v ?? rec.flat) as unknown;
    if (nr === 4 && nc === 4 && Array.isArray(vals) && vals.length === 16) {
      const v = (vals as unknown[]).map(n => clamp01(Number(n)));
      return [v.slice(0,4), v.slice(4,8), v.slice(8,12), v.slice(12,16)];
    }
  }
  // flat [16]
  if (Array.isArray(md) && md.length === 16) {
    const v = (md as unknown[]).map(n => clamp01(Number(n)));
    return [v.slice(0,4), v.slice(4,8), v.slice(8,12), v.slice(12,16)];
  }
  return null;
}

/** Normalize a map of matrices into Record<string, number[][]>. */
function normalizeTransMap(trans: Record<string, unknown> | null): Record<string, number[][]> {
  const out: Record<string, number[][]> = {};
  if (!trans || !isRec(trans)) return out;
  for (const [k, v] of Object.entries(trans)) {
    const m = extractRows4x4(v);
    if (m) out[k] = m;
  }
  return out;
}

function makeHistTrace(x: number[], name: string): Partial<PlotData> {
  return {
    type: "histogram",
    x,
    name,
    marker: { line: { width: 1 } },
    xbins: XBINS_0_1, // force [0,1] so single values render predictably
    hovertemplate: "p=%{x:.3f}<br>count=%{y}<extra></extra>",
  };
}

/* --------- component --------- */

export function ProbHists({
  root,
  trans,
  rootKey, // <-- pass the API's root_name here to exclude it
}: {
  root: number[] | null;
  trans: Record<string, unknown> | null;
  rootKey?: string | null;
}) {
  // sanitize root
  const rootClean = useMemo(
    () => (Array.isArray(root) ? root.map(n => clamp01(Number(n))).filter(n => Number.isFinite(n)) : []),
    [root]
  );

  // normalize and optionally exclude the root matrix
  const transMap = useMemo(() => {
    const m = normalizeTransMap(trans);
    if (rootKey && rootKey in m) {
      const { [rootKey]: _omit, ...rest } = m;
      return rest;
    }
    return m;
  }, [trans, rootKey]);

  // --- NEW: selector for a single entry aggregated across matrices ---
  const [rowIdx, setRowIdx] = useState<number>(2); // default: "array 2"
  const [colIdx, setColIdx] = useState<number>(3); // default: "element 3"

  const selectedCellValues = useMemo(() => {
    const r = Math.max(1, Math.min(4, rowIdx)) - 1;
    const c = Math.max(1, Math.min(4, colIdx)) - 1;
    const vals: number[] = [];
    for (const m of Object.values(transMap)) {
      const v = m[r]?.[c];
      if (Number.isFinite(v)) vals.push(clamp01(Number(v)));
    }
    return vals;
  }, [transMap, rowIdx, colIdx]);

  // Also keep your 16-per-cell grid histogram collection if you still want it:
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

  const rootOK = useMemo(() => rootClean.length === 4, [rootClean]);
  const anyTrans = useMemo(() => Object.keys(transMap).length > 0, [transMap]);

  return (
    <div className="space-y-6">
      {/* Root 4 histograms */}
      <section>
        <h2 className="text-lg font-semibold mb-2">Root probability (final)</h2>
        {rootOK ? (
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            {rootClean.map((val, i) => (
              <div key={`root-${i}`} className="border rounded p-2">
                <Plot
                  data={[makeHistTrace([val], `π${i + 1}`)]}
                  layout={{ ...smallLayout, title: { text: `π${i + 1}` }, bargap: 0.05 }}
                  config={baseConfig}
                  style={{ width: "100%", height: "100%" }}
                />
                <div className="text-xs text-gray-500 mt-1">
                  value: <span className="font-mono">{val.toFixed(6)}</span>
                </div>
              </div>
            ))}
          </div>
        ) : (
          <div className="text-sm text-gray-600">No root_prob_final (expected 4 values).</div>
        )}
      </section>

      {/* NEW: One-cell distribution across matrices (excluding root) */}
      <section>
        <h2 className="text-lg font-semibold mb-2">Transition entry distribution (across nodes)</h2>
        {!anyTrans ? (
          <div className="text-sm text-gray-600">No transition matrices available.</div>
        ) : (
          <div className="border rounded p-3 space-y-3">
            <div className="flex flex-wrap items-center gap-3">
              <label className="text-sm">
                Row&nbsp;
                <select
                  className="border rounded px-2 py-1 text-sm"
                  value={rowIdx}
                  onChange={(e) => setRowIdx(Number(e.target.value))}
                >
                  {[1, 2, 3, 4].map(n => <option key={n} value={n}>{n}</option>)}
                </select>
              </label>
              <label className="text-sm">
                Col&nbsp;
                <select
                  className="border rounded px-2 py-1 text-sm"
                  value={colIdx}
                  onChange={(e) => setColIdx(Number(e.target.value))}
                >
                  {[1, 2, 3, 4].map(n => <option key={n} value={n}>{n}</option>)}
                </select>
              </label>
              {rootKey ? (
                <div className="text-xs text-gray-600">Excluding root matrix: <code>{rootKey}</code></div>
              ) : null}
              <div className="text-xs text-gray-600">n = {selectedCellValues.length}</div>
            </div>

            <Plot
              data={[makeHistTrace(selectedCellValues, `m${rowIdx}${colIdx}`)]}
              layout={{ ...bigLayout, title: { text: `m${rowIdx}${colIdx} across nodes (n=${selectedCellValues.length})` }, bargap: 0.05 }}
              config={baseConfig}
              style={{ width: "100%", height: "100%" }}
            />
          </div>
        )}
      </section>

      {/* Transition 16 histograms (per-cell small multiples) */}
      <section>
        <h2 className="text-lg font-semibold mb-2">Transition matrix entries (final)</h2>
        {transCells.some(arr => arr.length) ? (
          <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-4 gap-4">
            {transCells.map((arr, idx) => {
              const i = Math.floor(idx / 4) + 1;
              const j = (idx % 4) + 1;
              return (
                <div key={`m-${idx}`} className="border rounded p-2">
                  <Plot
                    data={[makeHistTrace(arr, `m${i}${j}`)]}
                    layout={{ ...smallLayout, title: { text: `m${i}${j} (n=${arr.length})` }, bargap: 0.05 }}
                    config={baseConfig}
                    style={{ width: "100%", height: "100%" }}
                  />
                </div>
              );
            })}
          </div>
        ) : (
          <div className="text-sm text-gray-600">
            No transition probabilities found in <code>trans_prob_final</code>.
          </div>
        )}
      </section>
    </div>
  );
}
