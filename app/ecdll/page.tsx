// app/ecdll/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useState } from "react";
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
 * - alpha in (0,1): compresses big early gaps, stretches late small gaps.
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
  const [dirSum, setDirSum] = useState<SummaryRow | null>(null);
  const [parSum, setParSum] = useState<SummaryRow | null>(null);
  const [sshSum, setSshSum] = useState<SummaryRow | null>(null);

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

  /* -------- Build traces -------- */

  const compositeTraces = useMemo(() => {
    const traces: Partial<PlotData>[] = [];
    const push = (name: string, pts: Pt[] | null) => {
      if (!pts || !pts.length) return;
      const t = transformGapPower(pts, alpha);
      traces.push({
        type: "scatter",
        mode: "lines",
        name,
        x: t.map((p) => p[0]),
        y: t.map((p) => p[1]),
      });
    };
    push("Dirichlet", dirPts);
    push("Parsimony", parPts);
    push("SSH", sshPts);
    return traces;
  }, [dirPts, parPts, sshPts, alpha]);

  const makeTrace = (name: string, pts: Pt[] | null): Partial<PlotData> | null =>
    pts && pts.length
      ? { type: "scatter", mode: "lines", name, x: pts.map((p) => p[0]), y: pts.map((p) => p[1]) }
      : null;

  const dirTrace = useMemo(() => makeTrace("Dirichlet", dirPts), [dirPts]);
  const parTrace = useMemo(() => makeTrace("Parsimony", parPts), [parPts]);
  const sshTrace = useMemo(() => makeTrace("SSH", sshPts), [sshPts]);

  return (
    <div className="p-6 max-w-[1400px] mx-auto space-y-5">
      <div className="flex items-end gap-3">
        <h1 className="text-2xl font-bold">ECD-LL Convergence</h1>
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
          className="ml-auto px-4 py-2 rounded bg-black text-white disabled:opacity-50"
          disabled={!job || rep === null || loading}
        >
          {loading ? "Loading…" : "Reload"}
        </button>
      </div>

      {/* Summary table (for the selected repetition) */}
      <section className="border rounded p-3">
        <h2 className="text-lg font-semibold mb-2">Summary (rep {rep ?? "–"})</h2>
        <div className="overflow-x-auto">
          <table className="min-w-full text-sm">
            <thead>
              <tr className="text-left border-b">
                <th className="py-2 pr-4">Method</th>
                <th className="py-2 pr-4">Root</th>
                <th className="py-2 pr-4">Num Iterations</th>
                <th className="py-2 pr-4">LL Init</th>
                <th className="py-2 pr-4">ECD-LL First</th>
                <th className="py-2 pr-4">ECD-LL Final</th>
                <th className="py-2 pr-4">LL Final</th>
              </tr>
            </thead>
            <tbody>
              {[
                ["Dirichlet", dirSum] as const,
                ["Parsimony", parSum] as const,
                ["SSH", sshSum] as const,
              ].map(([name, s]) => (
                <tr key={name} className="border-b">
                  <td className="py-2 pr-4">{name}</td>
                  <td className="py-2 pr-4">{s?.root ?? "—"}</td>
                  <td className="py-2 pr-4">{s?.num_iterations ?? "—"}</td>
                  <td className="py-2 pr-4">
                    {s?.ll_init != null ? s.ll_init.toFixed(3) : "—"}
                  </td>
                  <td className="py-2 pr-4">
                    {s?.ecd_ll_first != null ? s.ecd_ll_first.toFixed(3) : "—"}
                  </td>
                  <td className="py-2 pr-4">
                    {s?.ecd_ll_final != null ? s.ecd_ll_final.toFixed(3) : "—"}
                  </td>
                  <td className="py-2 pr-4">
                    {s?.ll_final != null ? s.ll_final.toFixed(3) : "—"}
                  </td>
                </tr>
              ))}
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
        {dirPts && dirPts.length ? (
          <Plot
            data={[
              {
                type: "scatter",
                mode: "lines",
                name: "Dirichlet",
                x: dirPts.map((p) => p[0]),
                y: dirPts.map((p) => p[1]),
              } as Partial<PlotData>,
            ]}
            layout={baseLayout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
          />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      <section className="border rounded p-3">
        <h2 className="text-lg font-semibold mb-2">Parsimony</h2>
        {parPts && parPts.length ? (
          <Plot
            data={[
              {
                type: "scatter",
                mode: "lines",
                name: "Parsimony",
                x: parPts.map((p) => p[0]),
                y: parPts.map((p) => p[1]),
              } as Partial<PlotData>,
            ]}
            layout={baseLayout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
          />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      <section className="border rounded p-3">
        <h2 className="text-lg font-semibold mb-2">SSH</h2>
        {sshPts && sshPts.length ? (
          <Plot
            data={[
              {
                type: "scatter",
                mode: "lines",
                name: "SSH",
                x: sshPts.map((p) => p[0]),
                y: sshPts.map((p) => p[1]),
              } as Partial<PlotData>,
            ]}
            layout={baseLayout}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: "100%", height: "100%" }}
          />
        ) : (
          <div className="text-sm text-gray-600">No ECD-LL data.</div>
        )}
      </section>

      {err ? <div className="p-3 border rounded text-red-600">{err}</div> : null}
    </div>
  );
}
