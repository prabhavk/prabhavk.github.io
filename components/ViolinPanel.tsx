// app/components/ViolinPanel.tsx
"use client";

import React, { useEffect, useMemo, useRef, useState } from "react";
import dynamic from "next/dynamic";
import type * as Plotly from "plotly.js";

// Dynamically import Plotly to avoid SSR issues and type the component
type PlotParams = import("react-plotly.js").PlotParams;
const Plot = dynamic(() => import("react-plotly.js"), { ssr: false }) as unknown as React.ComponentType<PlotParams>;

const METHOD_ORDER = ["Parsimony", "Dirichlet", "SSH"] as const;
type Method = (typeof METHOD_ORDER)[number];

type Row = { method: Method; ll_final: number };
type ApiResp =
  | { data: Array<{ method: string; ll_final: number }> }
  | { error: string };

function normalizeMethod(m: unknown): Method | null {
  const s = String(m ?? "").trim().toUpperCase();
  if (s === "PARSIMONY") return "Parsimony";
  if (s === "DIRICHLET") return "Dirichlet";
  if (s === "SSH") return "SSH";
  return null;
}

function isAbortError(e: unknown): e is DOMException {
  return e instanceof DOMException && e.name === "AbortError";
}

export default function ViolinPanel({
  initialRoot,
  allRoots,
}: {
  initialRoot: string;
  allRoots: string[];
}) {
  const [root, setRoot] = useState(initialRoot);
  const [rows, setRows] = useState<Row[] | null>(null);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);
  const abortRef = useRef<AbortController | null>(null);

  async function fetchData(r: string) {
    // cancel any in-flight request
    abortRef.current?.abort();
    const ctrl = new AbortController();
    abortRef.current = ctrl;

    setLoading(true);
    setErr(null);
    try {
      const res = await fetch(`/api/llfinal?root=${encodeURIComponent(r)}`, {
        cache: "no-store",
        signal: ctrl.signal,
      });
      const json: ApiResp = await res.json();

      if (!res.ok || "error" in json) {
        throw new Error(
          "error" in json ? json.error : `HTTP ${res.status} while loading data`
        );
      }

      const normalized: Row[] = (json.data ?? [])
        .map((d) => {
          const m = normalizeMethod(d.method);
          const y = Number(d.ll_final);
          return m && Number.isFinite(y) ? { method: m, ll_final: y } : null;
        })
        .filter((x): x is Row => x !== null);

      setRows(normalized);
    } catch (e: unknown) {
      if (isAbortError(e)) return; // ignore aborted fetches
      const msg = e instanceof Error ? e.message : "Failed to load data";
      setErr(msg);
      setRows([]);
    } finally {
      setLoading(false);
    }
  }

  useEffect(() => {
    fetchData(root);
    // cleanup on unmount
    return () => abortRef.current?.abort();
  }, [root]);

  const traces: Partial<Plotly.Data>[] = useMemo(() => {
    // Always return 3 traces, even if some are empty, so x-axis order stays fixed
    const byMethod: Record<Method, number[]> = {
      Parsimony: [],
      Dirichlet: [],
      SSH: [],
    };
    for (const row of rows ?? []) {
      byMethod[row.method].push(row.ll_final);
    }

    return METHOD_ORDER.map((m): Partial<Plotly.Data> => ({
      type: "violin",
      name: m,
      y: byMethod[m],
      box: { visible: true },
      meanline: { visible: true },
      points: "all",
      hovertemplate: "<b>%{fullData.name}</b><br>ll_final=%{y}<extra></extra>",
    }));
  }, [rows]);

  const layout: Partial<Plotly.Layout> = useMemo(
    () => ({
      title: { text: `ll_final distributions @ root = ${root}` },
      xaxis: {
        title: { text: "Initialization method" },
        categoryorder: "array",
        // Plotly typings expect string[], so widen from readonly tuple:
        categoryarray: METHOD_ORDER as unknown as string[],
      },
      yaxis: { title: { text: "ll_final" } },
      showlegend: true,
      margin: { l: 50, r: 20, t: 50, b: 50 },
      violinmode: "group", // side-by-side violins
    }),
    [root]
  );

  return (
    <div className="space-y-4">
      <div className="flex items-center gap-3">
        <label className="text-sm font-medium">Root node</label>
        <select
          className="border rounded-lg px-2 py-1"
          value={root}
          onChange={(e) => setRoot(e.target.value)}
        >
          {allRoots.map((r) => (
            <option key={r} value={r}>
              {r}
            </option>
          ))}
        </select>
        {loading && <span className="text-sm opacity-70">Loading…</span>}
        {err && <span className="text-sm text-red-600">{err}</span>}
      </div>

      <Plot
        data={traces}
        layout={layout}
        config={{ responsive: true, displaylogo: false }}
        style={{ width: "100%", height: "500px" }}
      />
    </div>
  );
}
