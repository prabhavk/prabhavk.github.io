// app/pattern-weights/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import dynamic from "next/dynamic";
import type { Data, Layout } from "plotly.js";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

type ApiRows = {
  job_id: string;
  method: string; // normalized: "dirichlet" | "parsimony" | "hss"
  n: number;
  site_index: number[];
  weight: number[];
  cum_value: number[];
};

type MethodKey = "dirichlet" | "parsimony" | "hss";
type MethodLabel = "Dirichlet" | "Parsimony" | "HSS";

const METHOD_LABELS: Record<MethodKey, MethodLabel> = {
  dirichlet: "Dirichlet",
  parsimony: "Parsimony",
  hss: "HSS",
};

const COLORS: Record<MethodKey, string> = {
  dirichlet: "#D2691E",
  parsimony: "#FF6B3D",
  hss: "#BB1E10",
};

function labelToKey(label: MethodLabel): MethodKey {
  switch (label) {
    case "Dirichlet": return "dirichlet";
    case "Parsimony": return "parsimony";
    case "HSS":       return "hss";
  }
  // exhaustive, but TS needs a fallback:
  return "dirichlet";
}

function formatInt(n: number): string {
  return new Intl.NumberFormat().format(n);
}

export default function PatternWeightsPage() {
  const [job, setJob] = useState<string>("");
  const [methodLabel, setMethodLabel] = useState<MethodLabel>("Dirichlet");
  const method: MethodKey = useMemo(() => labelToKey(methodLabel), [methodLabel]);

  const [data, setData] = useState<ApiRows | null>(null);
  const [loading, setLoading] = useState<boolean>(false);
  const [err, setErr] = useState<string | null>(null);

  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId") ?? "";
      setJob(saved);
    } catch {/* ignore */}
  }, []);

  const fetchRows = useCallback(async (jid: string, mk: MethodKey) => {
    const u = new URL("/api/pattern-weights", window.location.origin);
    u.searchParams.set("job", jid);
    u.searchParams.set("method", mk);
    const res = await fetch(u.toString(), { cache: "no-store" });
    const j = (await res.json()) as unknown;
    if (!res.ok) {
      const msg = (j as { error?: string })?.error ?? `HTTP ${res.status}`;
      throw new Error(msg);
    }
    const rows = j as ApiRows;
    if (!Array.isArray(rows.site_index) || !Array.isArray(rows.weight) || !Array.isArray(rows.cum_value)) {
      throw new Error("Malformed response");
    }
    return rows;
  }, []);

  const load = useCallback(async () => {
    if (!job) {
      setData(null);
      setErr("No job selected");
      return;
    }
    setLoading(true);
    setErr(null);
    try {
      const rows = await fetchRows(job, method);
      setData(rows);
    } catch (e) {
      const msg = e instanceof Error ? e.message : "Failed to load";
      setErr(msg);
      setData(null);
    } finally {
      setLoading(false);
    }
  }, [fetchRows, job, method]);

  useEffect(() => { void load(); }, [load]);

  const traces: Partial<Data>[] = useMemo(() => {
    if (!data || !data.n) return [];
    const x = data.site_index;
    const yW = data.weight;
    const yC = data.cum_value;
    const col = COLORS[method];

    const tWeight: Partial<Data> = {
      type: "scatter",
      mode: "lines",
      name: "weight",
      x,
      y: yW,
      line: { color: col, width: 2 },
      hovertemplate: "site=%{x}<br>weight=%{y:.6f}<extra></extra>",
    };

    const tCum: Partial<Data> = {
      type: "scatter",
      mode: "lines",
      name: "cum_value",
      x,
      y: yC,
      line: { color: "black", width: 2, dash: "dot" },
      hovertemplate: "site=%{x}<br>cum=%{y:.6f}<extra></extra>",
      yaxis: "y2",
    };

    return [tWeight, tCum];
  }, [data, method]);

  const layout: Partial<Layout> = useMemo(() => {
    const titleText = data
      ? `Pattern Weights · ${METHOD_LABELS[method]} · job=${data.job_id}`
      : `Pattern Weights · ${METHOD_LABELS[method]}`;

    return {
      title: { text: titleText },
      height: 560,
      margin: { l: 60, r: 60, t: 60, b: 60 },
      legend: { orientation: "h" },
      xaxis: {
        title: { text: "site index" },     // ← must be an object
        zerolinecolor: "#e5e7eb",
        gridcolor: "#f3f4f6",
      },
      yaxis: {
        title: { text: "weight" },         // ← must be an object
        zerolinecolor: "#e5e7eb",
        gridcolor: "#f3f4f6",
      },
      yaxis2: {
        title: { text: "cum_value" },      // ← must be an object
        overlaying: "y",
        side: "right",
        zerolinecolor: "#e5e7eb",
        gridcolor: "#f3f4f6",
      },
      paper_bgcolor: "white",
      plot_bgcolor: "white",
    };
  }, [data, method]);

  const exportCsv = useCallback(() => {
    if (!data || !data.n) return;
    const lines = ["site_index,weight,cum_value"];
    for (let i = 0; i < data.site_index.length; i++) {
      lines.push([data.site_index[i], data.weight[i], data.cum_value[i]].join(","));
    }
    const content = lines.join("\n");
    const blob = new Blob([content], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const fname = `${data.job_id}__pattern-weights__${data.method}.csv`;
    const a = document.createElement("a");
    a.href = url;
    a.download = fname;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }, [data]);

  return (
    <div className="p-6 max-w-[1200px] mx-auto">
      <h1 className="text-2xl font-bold mb-4">Pattern Weights</h1>

      <div className="flex flex-wrap gap-3 items-end mb-4">
        <label className="flex flex-col">
          <span className="text-sm">Job ID</span>
          <input
            className="border rounded px-3 py-2 min-w-[340px]"
            placeholder="job-…"
            value={job}
            onChange={(e) => setJob(e.target.value)}
          />
        </label>

        <label className="flex flex-col">
          <span className="text-sm">Method</span>
          <select
            className="border rounded px-3 py-2"
            value={methodLabel}
            onChange={(e) => setMethodLabel(e.target.value as MethodLabel)}
          >
            {(["Dirichlet", "Parsimony", "HSS"] as MethodLabel[]).map((m) => (
              <option key={m} value={m}>{m}</option>
            ))}
          </select>
        </label>

        <button
          className="px-3 h-9 rounded bg-black text-white disabled:opacity-50"
          onClick={load}
          disabled={!job || loading}
          title="Reload"
        >
          {loading ? "Loading…" : "Reload"}
        </button>

        <div className="ml-auto text-sm">
          {err ? <span className="text-red-600">{err}</span> : null}
          {data && !err ? (
            <span className="text-gray-700">
              {formatInt(data.n)} rows · method <span className="font-mono">{data.method}</span>
            </span>
          ) : null}
        </div>
      </div>

      <section className="bg-white border rounded p-3">
        {!data || !data.n ? (
          <div className="p-3 border rounded">No data</div>
        ) : (
          <>
            <div className="flex items-center justify-between mb-2">
              <h2 className="text-lg font-semibold">Weights &amp; Cumulative</h2>
              <button
                className="px-3 py-1 text-sm rounded bg-black text-white disabled:opacity-50"
                onClick={exportCsv}
                disabled={!data || !data.n}
              >
                Download CSV
              </button>
            </div>
            <Plot
              data={traces}
              layout={layout}
              config={{ displayModeBar: false, responsive: true }}
              style={{ width: "100%", height: "100%" }}
            />

            <div className="mt-4 text-sm text-gray-700">
              <p>
                <span className="font-semibold">Tip:</span> This plot draws per-site{" "}
                <span className="font-mono">weight</span> (left axis) and{" "}
                <span className="font-mono">cum_value</span> (right axis) for{" "}
                <span className="font-mono">{METHOD_LABELS[method]}</span>.
              </p>
            </div>
          </>
        )}
      </section>
    </div>
  );
}
