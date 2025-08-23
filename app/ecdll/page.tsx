// app/ecdll/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
import dynamic from "next/dynamic";
import type { Layout, PlotData } from "plotly.js";
import type { PlotParams } from "react-plotly.js";

const Plot = dynamic<PlotParams>(() => import("react-plotly.js"), { ssr: false });

/* ---------------- Types ---------------- */

type Pt = [number, number]; // [iter, ecd_ll]

type ApiOk = {
  ok: true;
  job_id: string;
  method: string;
  points: Pt[];
};
type ApiErr = { ok: false; error: string };

type MethodKey = "dirichlet" | "parsimony" | "ssh";
const METHOD_LABEL: Record<MethodKey, string> = {
  dirichlet: "Dirichlet",
  parsimony: "Parsimony",
  ssh: "SSH",
};

type SummaryRow = {
  ok: true;
  method: MethodKey;
  root: number[] | null;
  num_iterations: number;
  initial_ll: number | null;
  final_ll: number | null;
  ecd_ll_first: number | null;
  ecd_ll_final: number | null;  
};
type SummaryResp =
  | { ok: true; job_id: string; rows: SummaryRow[] }
  | { ok: false; error: string };

/* ---------------- Helpers ---------------- */

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

/** Gap^alpha transform for composite plot */
function transformGapPower(points: Pt[], alpha: number = 0.5, eps = 1e-9): Pt[] {
  if (!Array.isArray(points) || points.length === 0) return points;
  const ys = points.map(([, ll]) => ll).filter((v) => Number.isFinite(v));
  if (ys.length === 0) return points;
  const best = Math.max(...ys);
  const gaps = points.map(([, ll]) => best - ll);
  const gMax = Math.max(...gaps.filter((g) => Number.isFinite(g)), 0);

  return points.map(([it, _ll], i) => {
    const g = gaps[i];
    if (!Number.isFinite(g)) return [it, 0];
    const y = Math.pow(g / (gMax + eps), alpha);
    return [it, y];
  });
}

function transformLogGap(points: Pt[], eps = 1e-12): Pt[] {
  if (!points.length) return points;
  const ys = points.map(([, ll]) => ll).filter((v) => Number.isFinite(v));
  if (!ys.length) return points;
  const best = Math.max(...ys);
  return points.map(([it, ll]) => [it, Math.log(best - ll + eps)]);
}

const STYLE: Record<MethodKey, Partial<PlotData>> = {
  dirichlet: { line: { dash: "solid", width: 2 } },
  parsimony: { line: { dash: "dash", width: 2 } },
  ssh: { line: { dash: "dot", width: 2 } },
};

/* ---------------- Page ---------------- */

export default function EcdllPage() {
  const [job, setJob] = useState<string>("");
  const [err, setErr] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  // Raw points per method
  const [dirPts, setDirPts] = useState<Pt[] | null>(null);
  const [parPts, setParPts] = useState<Pt[] | null>(null);
  const [sshPts, setSshPts] = useState<Pt[] | null>(null);

  // Summary rows
  const [summary, setSummary] = useState<Record<MethodKey, SummaryRow | null>>({
    dirichlet: null,
    parsimony: null,
    ssh: null,
  });

  // Tuning for composite transform
  const [alpha, setAlpha] = useState<number>(0.5);

  // Read selected job id
  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId");
      if (saved) setJob(saved);
    } catch {
      /* ignore */
    }
  }, []);

  const fetchOne = useCallback(
    async (method: MethodKey): Promise<Pt[]> => {
      const u = new URL("/api/ecdll", window.location.origin);
      u.searchParams.set("job_id", job);
      u.searchParams.set("method", method);

      const res = await fetch(u.toString(), { cache: "no-store" });
      const ct = res.headers.get("content-type") || "";
      if (!ct.toLowerCase().includes("application/json")) {
        const body = await res.text().catch(() => "");
        throw new Error(
          `Expected JSON for ${method}, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`
        );
      }
      const j = (await res.json()) as ApiOk | ApiErr;

      if (!res.ok || !("ok" in j) || (j as ApiErr).ok === false) {
        const msg = "ok" in j && (j as ApiErr).ok === false ? (j as ApiErr).error : `HTTP ${res.status}`;
        throw new Error(`${METHOD_LABEL[method]}: ${msg}`);
      }
      return normalizePoints((j as ApiOk).points);
    },
    [job]
  );

  const fetchSummary = useCallback(async (): Promise<Record<MethodKey, SummaryRow | null>> => {
    const u = new URL("/api/ecdll/summary", window.location.origin);
    u.searchParams.set("job_id", job);

    const res = await fetch(u.toString(), { cache: "no-store" });
    const ct = res.headers.get("content-type") || "";
    if (!ct.toLowerCase().includes("application/json")) {
      const body = await res.text().catch(() => "");
      throw new Error(`Expected JSON for summary, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`);
    }
    const j = (await res.json()) as SummaryResp;
    if (!res.ok || j.ok === false) {
      const msg = j.ok === false ? j.error : `HTTP ${res.status}`;
      throw new Error(msg);
    }

    const map: Record<MethodKey, SummaryRow | null> = {
      dirichlet: null,
      parsimony: null,
      ssh: null,
    };
    for (const r of j.rows) {
      map[r.method] = r;
    }
    return map;
  }, [job]);

  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Choose a job on the Precomputed Results page.");
      setDirPts(null);
      setParPts(null);
      setSshPts(null);
      setSummary({ dirichlet: null, parsimony: null, ssh: null });
      return;
    }
    setErr(null);
    setLoading(true);
    try {
      const [d, p, s, sum] = await Promise.allSettled([
        fetchOne("dirichlet"),
        fetchOne("parsimony"),
        fetchOne("ssh"),
        fetchSummary(),
      ]);

      setDirPts(d.status === "fulfilled" ? d.value : []);
      setParPts(p.status === "fulfilled" ? p.value : []);
      setSshPts(s.status === "fulfilled" ? s.value : []);
      setSummary(sum.status === "fulfilled" ? sum.value : { dirichlet: null, parsimony: null, ssh: null });

      const errs: string[] = [];
      if (d.status === "rejected") errs.push(d.reason?.message || "Dirichlet fetch failed");
      if (p.status === "rejected") errs.push(p.reason?.message || "Parsimony fetch failed");
      if (s.status === "rejected") errs.push(s.reason?.message || "SSH fetch failed");
      if (sum.status === "rejected") errs.push(sum.reason?.message || "Summary fetch failed");
      if (errs.length) setErr(errs.join(" | "));
    } catch (e) {
      setErr(e instanceof Error ? e.message : "Failed to load ECD-LL data");
      setDirPts([]);
      setParPts([]);
      setSshPts([]);
      setSummary({ dirichlet: null, parsimony: null, ssh: null });
    } finally {
      setLoading(false);
    }
  }, [job, fetchOne, fetchSummary]);

  useEffect(() => {
    if (job) void load();
  }, [job, load]);

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

  /* -------- Build traces -------- */

  const compositeTraces = useMemo(() => {
    const traces: Partial<PlotData>[] = [];
    if (dirPts && dirPts.length) {
      const t = transformGapPower(dirPts, alpha);
      traces.push({ type: "scatter", mode: "lines", name: "Dirichlet", x: t.map((p) => p[0]), y: t.map((p) => p[1]), ...STYLE.dirichlet });
    }
    if (parPts && parPts.length) {
      const t = transformGapPower(parPts, alpha);
      traces.push({ type: "scatter", mode: "lines", name: "Parsimony", x: t.map((p) => p[0]), y: t.map((p) => p[1]), ...STYLE.parsimony });
    }
    if (sshPts && sshPts.length) {
      const t = transformGapPower(sshPts, alpha);
      traces.push({ type: "scatter", mode: "lines", name: "SSH", x: t.map((p) => p[0]), y: t.map((p) => p[1]), ...STYLE.ssh });
    }
    return traces;
  }, [dirPts, parPts, sshPts, alpha]);

  const dirTrace = useMemo<Partial<PlotData> | null>(() => {
    if (!dirPts || !dirPts.length) return null;
    return { type: "scatter", mode: "lines", name: "Dirichlet", x: dirPts.map((p) => p[0]), y: dirPts.map((p) => p[1]), ...STYLE.dirichlet };
  }, [dirPts]);

  const parTrace = useMemo<Partial<PlotData> | null>(() => {
    if (!parPts || !parPts.length) return null;
    return { type: "scatter", mode: "lines", name: "Parsimony", x: parPts.map((p) => p[0]), y: parPts.map((p) => p[1]), ...STYLE.parsimony };
  }, [parPts]);

  const sshTrace = useMemo<Partial<PlotData> | null>(() => {
    if (!sshPts || !sshPts.length) return null;
    return { type: "scatter", mode: "lines", name: "SSH", x: sshPts.map((p) => p[0]), y: sshPts.map((p) => p[1]), ...STYLE.ssh };
  }, [sshPts]);

  /* -------- Table helpers -------- */

  const formatRoot = (root: number[] | null) =>
    root ? root.map((v) => v.toFixed(6)).join(", ") : "—";

  const S: MethodKey[] = ["dirichlet", "parsimony", "ssh"];

  return (
    <div className="p-6 max-w-[1400px] mx-auto space-y-5">
      <div className="flex items-end gap-3">
        <h1 className="text-2xl font-bold">ECD-LL Convergence</h1>
        <div className="text-sm ml-2">
          Job: <span className="font-mono">{job || "(none)"} </span>
        </div>
        <button
          type="button"
          onClick={() => void load()}
          className="ml-auto px-4 py-2 rounded bg-black text-white disabled:opacity-50"
          disabled={!job || loading}
        >
          {loading ? "Loading…" : "Reload"}
        </button>
      </div>

      {/* Summary table */}
      <section className="border rounded p-3">
        <h2 className="text-lg font-semibold mb-3">Run summary</h2>
        <div className="overflow-auto">
          <table className="min-w-full text-sm">
            <thead className="text-left">
              <tr className="border-b">
                <th className="py-2 pr-3">Method</th>
                <th className="py-2 pr-3">root</th>
                <th className="py-2 pr-3">num_iterations</th>
                <th className="py-2 pr-3">initial_ll</th>
                <th className="py-2 pr-3">final_ll</th>
                <th className="py-2 pr-3">ecd_ll_first</th>
                <th className="py-2 pr-3">ecd_ll_final</th>                
              </tr>
            </thead>
            <tbody>
              {S.map((m) => {
                const r = summary[m];
                return (
                  <tr key={m} className="border-b last:border-0">
                    <td className="py-2 pr-3 font-medium">{METHOD_LABEL[m]}</td>
                    <td className="py-2 pr-3 font-mono">{r ? formatRoot(r.root) : "—"}</td>
                    <td className="py-2 pr-3">{r?.num_iterations ?? "—"}</td>
                    <td className="py-2 pr-3 font-mono">{r?.initial_ll ?? "—"}</td>
                    <td className="py-2 pr-3 font-mono">{r?.final_ll ?? "—"}</td>
                    <td className="py-2 pr-3 font-mono">{r?.ecd_ll_first ?? "—"}</td>
                    <td className="py-2 pr-3 font-mono">{r?.ecd_ll_final ?? "—"}</td>                    
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
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
        <h2 className="text-lg font-semibold mb-2">Composite (Gap^α transform)</h2>
        {compositeTraces.length ? (
          <Plot
            data={compositeTraces}
            layout={compositeLayout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
          />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data found for any method.</div>
        )}
      </section>

      {/* Individual plots */}
      <section className="border rounded p-3">
        <h2 className="text-lg font-semibold mb-2">Dirichlet</h2>
        {dirTrace ? (
          <Plot data={[dirTrace]} layout={baseLayout} config={{ displayModeBar: false, responsive: true }} style={{ width: "100%", height: "100%" }} />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      <section className="border rounded p-3">
        <h2 className="text-lg font-semibold mb-2">Parsimony</h2>
        {parTrace ? (
          <Plot data={[parTrace]} layout={baseLayout} config={{ displayModeBar: false, responsive: true }} style={{ width: "100%", height: "100%" }} />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      <section className="border rounded p-3">
        <h2 className="text-lg font-semibold mb-2">SSH</h2>
        {sshTrace ? (
          <Plot data={[sshTrace]} layout={baseLayout} config={{ displayModeBar: false, responsive: true }} style={{ width: "100%", height: "100%" }} />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      {err ? <div className="p-3 border rounded text-red-600">{err}</div> : null}
    </div>
  );
}
