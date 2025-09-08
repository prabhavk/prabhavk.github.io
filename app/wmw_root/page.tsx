// app/wmw_comp_root/page.tsx
"use client";

const showNode = (s: string) => s.replace(/^h_/, "");

import React, { useEffect, useMemo, useState } from "react";

/** ---- Types that match /api/wmw_comp_root (with normalization) ---- */
type MethodName = "Parsimony" | "Dirichlet" | "HSS";

type RawApiResp = {
  job_id: string;
  method: string;
  alpha: number;
  nodes: string[];
  sizes: Record<string, number>;
  p_value: number[][];
  z: number[][];
  nA: number[][];
  nB: number[][];
  metric?: "final" | "init";
};

type ApiResp = Omit<RawApiResp, "method"> & { method: MethodName };

/** ---------- robust fetch helper ---------- */
function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isApiErr(x: unknown): x is { error: string } {
  return isRecord(x) && typeof x["error"] === "string";
}
async function fetchJSON<T>(url: string, init?: RequestInit): Promise<T> {
  const res = await fetch(url, { ...init, headers: { accept: "application/json", ...(init?.headers || {}) } });
  const ct = (res.headers.get("content-type") || "").toLowerCase();
  const body = await res.text();

  if (!ct.includes("application/json")) {
    throw new Error(`Expected JSON but got ${ct || "unknown"} (HTTP ${res.status}). ` + body.slice(0, 200));
  }
  let json: unknown;
  try {
    json = body ? JSON.parse(body) : {};
  } catch (e) {
    throw new Error(`Invalid JSON: ${(e as Error).message}`);
  }
  if (!res.ok) {
    const msg = isApiErr(json) ? json.error : `HTTP ${res.status}`;
    throw new Error(msg);
  }
  return json as T;
}

/** Normalize method labels coming from server */
function normMethod(s: string): MethodName | null {
  const m = s.trim().toLowerCase();
  if (m.includes("pars")) return "Parsimony";
  if (m.startsWith("dir") || m.includes("dirich")) return "Dirichlet";
  if (m.includes("hss") || m === "ssh") return "HSS";
  return null;
}

/* ---------------------------------------- */

const ALL_METHODS: MethodName[] = ["Parsimony", "Dirichlet", "HSS"];
type Metric = "final" | "init";

type Bundle = {
  data: ApiResp;
  nodes: string[];
  sizes: Record<string, number>;
  agreeCount: number[][];
  available: MethodName[];
};

/** ---------- CSV helpers ---------- */
function csvEscape(value: unknown): string {
  const s = value == null ? "" : String(value);
  if (/[",\n\r\t]/.test(s)) {
    return `"${s.replace(/"/g, '""')}"`;
  }
  return s;
}
function rowsToCsv(rows: (string | number)[][]): string {
  return rows.map(r => r.map(csvEscape).join(",")).join("\n");
}
function downloadCsv(filename: string, rows: (string | number)[][]) {
  const blob = new Blob([rowsToCsv(rows)], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

/** ---------- Page ---------- */

export default function WmwCompRootPage() {
  const [jobId, setJobId] = useState<string>("");
  const [alpha, setAlpha] = useState<number>(0.05);
  const [method, setMethod] = useState<MethodName>("Parsimony");

  // metric bundles
  const [finalB, setFinalB] = useState<Bundle | null>(null);
  const [initB, setInitB] = useState<Bundle | null>(null);

  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);
  const [showZ, setShowZ] = useState<boolean>(false);

  // hydrate jobId
  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId");
      if (saved) setJobId(saved);
    } catch { /* ignore */ }
  }, []);

  const alphaStr = useMemo(() => alpha.toString(), [alpha]);

  async function fetchBundle(metric: Metric): Promise<Bundle> {
    // build URLs per requested method label
    const urls = ALL_METHODS.map((m) => {
      const u = new URL("/api/wmw_comp_root", window.location.origin);
      u.searchParams.set("job", jobId);
      u.searchParams.set("method", m);
      u.searchParams.set("alpha", String(alpha));
      u.searchParams.set("metric", metric);
      return u.toString();
    });

    // be tolerant of partial failure
    const settled = await Promise.allSettled(urls.map((u) => fetchJSON<RawApiResp>(u, { cache: "no-store" })));

    // normalize & index by normalized method
    const byMethod = new Map<MethodName, ApiResp>();
    for (const s of settled) {
      if (s.status !== "fulfilled") continue;
      const nm = normMethod(s.value.method);
      if (!nm) continue;
      byMethod.set(nm, { ...s.value, method: nm });
    }

    const available = Array.from(byMethod.keys());
    if (available.length === 0) {
      throw new Error("No method data returned by the server for this job/α/metric.");
    }

    // choose base for node order:
    // prefer the user-selected method if available, else the first available
    const baseMethod: MethodName = available.includes(method) ? method : available[0];
    const base = byMethod.get(baseMethod)!;

    // intersection of node sets across available methods, preserving base order
    const sets = available.map((m) => new Set(byMethod.get(m)!.nodes));
    const nodes = base.nodes.filter((n) => sets.every((S) => S.has(n)));
    const N = nodes.length;

    // build per-method node -> index maps
    const idx: Record<MethodName, Record<string, number>> = {} as any;
    for (const m of available) {
      const resp = byMethod.get(m)!;
      idx[m] = Object.fromEntries(resp.nodes.map((n, i) => [n, i]));
    }

    // agreement counts across available methods (p<alpha & z>0 means column>row)
    const agreeCount: number[][] = Array.from({ length: N }, () => Array(N).fill(0));
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        if (i === j) continue;
        const rowName = nodes[i], colName = nodes[j];
        let agree = 0;
        for (const m of available) {
          const resp = byMethod.get(m)!;
          const im = idx[m][rowName];
          const jm = idx[m][colName];
          const p = resp.p_value[im]?.[jm];
          const z = resp.z[im]?.[jm];
          if (Number.isFinite(p) && Number.isFinite(z) && (p as number) < alpha && (z as number) > 0) {
            agree++;
          }
        }
        agreeCount[i][j] = agree;
      }
    }

    // sizes from the selected base method
    const sizes = byMethod.get(baseMethod)!.sizes;

    // selected data (if user's chosen method not available, return baseMethod’s matrix)
    const selected = byMethod.get(baseMethod)!;

    return { data: selected, nodes, sizes, agreeCount, available };
  }

  async function run() {
    if (!jobId) {
      setErr("No job selected. Pick a job on the Precomputed Results page.");
      return;
    }
    setLoading(true);
    setErr(null);
    setFinalB(null);
    setInitB(null);

    try {
      const [bFinal, bInit] = await Promise.all([fetchBundle("final"), fetchBundle("init")]);
      setFinalB(bFinal);
      setInitB(bInit);
    } catch (e: unknown) {
      setErr(e instanceof Error ? e.message : "Failed to fetch results");
    } finally {
      setLoading(false);
    }
  }

  // 🔁 Auto-run whenever jobId, method, or alpha change
  useEffect(() => {
    if (!jobId) return;
    void run();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [jobId, method, alpha]);

  // helper to show which methods were available
  const availNote = (b: Bundle | null) =>
    b ? `Available methods: ${b.available.join(", ")}` : "";

  return (
    <div className="p-6 max-w-[95vw] mx-auto space-y-6">
      <header className="space-y-1">
        <h1 className="text-2xl text-black font-bold">WMW (one-sided) — Node vs Node</h1>
        <p className="text-sm text-gray-600">
          Cell p-values are for H₁: <em>column</em> &gt; <em>row</em> (one-tailed).
          Agreement counts how many methods support this (p &lt; α and z &gt; 0).
          We show results for both <strong>final</strong> and <strong>initial</strong> log-likelihoods.
        </p>
      </header>

      <section className="flex flex-wrap items-end gap-3">
        <label className="flex flex-col">
          <span className="text-sm text-black">Method</span>
          <select
            value={method}
            onChange={(e) => setMethod(e.target.value as MethodName)}
            className="border rounded-md text-blue-700 px-3 py-2"
            disabled={loading}
          >
            {ALL_METHODS.map((m) => <option key={m} value={m}>{m}</option>)}
          </select>
        </label>

        <label className="flex flex-col">
          <span className="text-sm">α (significance)</span>
          <input
            type="number"
            step="0.001"
            min={0.001}
            max={0.5}
            value={alphaStr}
            onChange={(e) => {
              const v = Number(e.target.value);
              if (Number.isFinite(v)) setAlpha(v);
            }}
            className="border rounded-md px-3 py-2 w-28"
            disabled={loading}
          />
        </label>

        <div className="flex items-center gap-2">
          <input
            id="showZ"
            type="checkbox"
            checked={showZ}
            onChange={(e) => setShowZ(e.target.checked)}
            className="h-4 w-4"
            disabled={loading}
          />
          <label htmlFor="showZ" className="text-sm">Show z-scores</label>
        </div>

        <div className="text-sm text-gray-600 ml-auto">
          Job: <span className="font-mono">{jobId || "(none)"}</span>
          {loading && <span className="ml-3 animate-pulse">Loading…</span>}
        </div>
      </section>

      {err && (
        <div className="p-3 rounded-md bg-red-50 text-red-700 border border-red-200">
          {err}
        </div>
      )}

      {/* -------- FINAL (ll_final) -------- */}
      {!err && finalB && (
        <>
          <h2 className="text-xl font-semibold text-black mt-2">Final log-likelihood (ll_final)</h2>
          <div className="text-xs text-gray-500 mb-1">{availNote(finalB)}</div>
          <Legend alpha={finalB.data.alpha} />
          <MatrixTable
            title={`P/Z Matrix — ${finalB.data.method}`}
            nodes={finalB.nodes}
            sizes={finalB.sizes}
            p={finalB.data.p_value}
            z={finalB.data.z}
            alpha={finalB.data.alpha}
            showZ={showZ}
            jobId={jobId}
            method={finalB.data.method}
            metric="final"
          />
          <AgreementLegend />
          <CountsMatrixTable
            title="Agreement Count (methods supporting column > row)"
            nodes={finalB.nodes}
            sizes={finalB.sizes}
            counts={finalB.agreeCount}
            maxCount={finalB.available.length}
            jobId={jobId}
            method={finalB.data.method}
            metric="final"
            alpha={finalB.data.alpha}
          />
        </>
      )}

      {/* -------- INITIAL (ll_init) -------- */}
      {!err && initB && (
        <>
          <h2 className="text-xl font-semibold text-black mt-6">Initial log-likelihood (ll_init)</h2>
          <div className="text-xs text-gray-500 mb-1">{availNote(initB)}</div>
          <Legend alpha={initB.data.alpha} />
          <MatrixTable
            title={`P/Z Matrix — ${initB.data.method}`}
            nodes={initB.nodes}
            sizes={initB.sizes}
            p={initB.data.p_value}
            z={initB.data.z}
            alpha={initB.data.alpha}
            showZ={showZ}
            jobId={jobId}
            method={initB.data.method}
            metric="init"
          />
          <AgreementLegend />
          <CountsMatrixTable
            title="Agreement Count (methods supporting column > row)"
            nodes={initB.nodes}
            sizes={initB.sizes}
            counts={initB.agreeCount}
            maxCount={initB.available.length}
            jobId={jobId}
            method={initB.data.method}
            metric="init"
            alpha={initB.data.alpha}
          />
        </>
      )}

      {!err && !loading && !finalB && !initB && (
        <p className="text-gray-500 text-sm">Pick a job, then method/α; results load automatically.</p>
      )}
    </div>
  );
}

/** ---------- P/Z Matrix table (with CSV download & footer sums) ---------- */

function MatrixTable({
  title,
  nodes,
  sizes,
  p,
  z,
  alpha,
  showZ,
  jobId,
  method,
  metric,
}: {
  title: string;
  nodes: string[];
  sizes: Record<string, number>;
  p: number[][];
  z: number[][];
  alpha: number;
  showZ: boolean;
  jobId: string;
  method: MethodName;
  metric: "final" | "init";
}) {
  function cellClass(i: number, j: number): string {
    if (i === j) return "bg-gray-100 text-gray-400";
    const pij = p[i]?.[j];
    const zij = z[i]?.[j];
    if (!Number.isFinite(pij) || !Number.isFinite(zij)) return "";
    if (pij < alpha && zij > 0) return "bg-emerald-100";
    if (zij < 0) return "bg-rose-50";
    return "bg-white";
  }

  const fmtP = (x: number) => {
    if (!Number.isFinite(x)) return "";
    if (x === 0) return "0";
    if (x < 1e-3) return x.toExponential(2);
    return x.toFixed(4);
  };
  const fmtZ = (x: number) => (Number.isFinite(x) ? x.toFixed(2) : "");

  // column sums of support (p<alpha && z>0), excluding diagonal
  const N = nodes.length;
  const colSupportSums = Array.from({ length: N }, (_, j) =>
    p.reduce((acc, _row, i) => {
      if (i === j) return acc;
      const pij = p[i]?.[j];
      const zij = z[i]?.[j];
      return Number.isFinite(pij) && Number.isFinite(zij) && (pij as number) < alpha && (zij as number) > 0
        ? acc + 1
        : acc;
    }, 0)
  );

  // CSV builders (strip h_ in headers and row labels)
  const buildCsvHeader = (): (string | number)[][] => [
    ["row/col", ...nodes.map(showNode)],
    ["sample size (per node)", ...nodes.map((n) => sizes[n] ?? 0)],
  ];

  const downloadPcsv = () => {
    const rows: (string | number)[][] = [...buildCsvHeader()];
    for (let i = 0; i < N; i++) {
      const r: (string | number)[] = [showNode(nodes[i])];
      for (let j = 0; j < N; j++) {
        r.push(i === j ? "—" : fmtP(p[i]?.[j]));
      }
      rows.push(r);
    }
    rows.push(["column sum (p<α & z>0)", ...colSupportSums]);
    const fname = `wmw_pvals_${metric}_${method}_${jobId}_alpha-${alpha}.csv`;
    downloadCsv(fname, rows);
  };

  const downloadZcsv = () => {
    const rows: (string | number)[][] = [...buildCsvHeader()];
    for (let i = 0; i < N; i++) {
      const r: (string | number)[] = [showNode(nodes[i])];
      for (let j = 0; j < N; j++) {
        r.push(i === j ? "—" : fmtZ(z[i]?.[j]));
      }
      rows.push(r);
    }
    rows.push(["column sum (p<α & z>0)", ...colSupportSums]);
    const fname = `wmw_zscores_${metric}_${method}_${jobId}_alpha-${alpha}.csv`;
    downloadCsv(fname, rows);
  };

  return (
    <section className="space-y-2">
      <div className="flex items-center justify-between">
        <h2 className="text-lg font-semibold">{title}</h2>
        <div className="flex gap-2">
          <button
            type="button"
            onClick={downloadPcsv}
            className="px-2 py-1 text-xs rounded bg-black text-white"
            title="Download p-values CSV"
          >
            CSV (p-values)
          </button>
          <button
            type="button"
            onClick={downloadZcsv}
            className="px-2 py-1 text-xs rounded bg-black text-white"
            title="Download z-scores CSV"
          >
            CSV (z-scores)
          </button>
        </div>
      </div>

      <div className="overflow-auto border rounded-xl">
        <table className="min-w-max text-sm">
          <thead className="sticky top-0 bg-white z-10">
            <tr>
              <Th className="sticky left-0 bg-white z-20">row ↓ / col →</Th>
              {nodes.map((n) => (
                <Th key={n} className="text-center">{showNode(n)}</Th>
              ))}
            </tr>
            <tr>
              <Th className="sticky left-0 bg-white z-20">sample size (per node)</Th>
              {nodes.map((n) => (
                <Th key={`${n}-n`} className="text-center">{sizes[n] ?? 0}</Th>
              ))}
            </tr>
          </thead>
          <tbody>
            {nodes.map((rowName, i) => (
              <tr key={rowName} className="hover:bg-gray-50">
                <Td mono className="sticky left-0 bg-white z-10">{showNode(rowName)}</Td>
                {nodes.map((_, j) => (
                  <Td key={`${i}-${j}`} center className={cellClass(i, j)}>
                    {i === j ? "—" : (
                      <div className="leading-tight">
                        <div className="font-mono">{fmtP(p[i]?.[j])}</div>
                        {showZ && <div className="text-xs text-gray-500">z={fmtZ(z[i]?.[j])}</div>}
                      </div>
                    )}
                  </Td>
                ))}
              </tr>
            ))}
          </tbody>

          {/* Footer: column sums of supports (p<α ∧ z>0), excl. diagonal */}
          <tfoot className="bg-gray-50">
            <tr className="font-semibold">
              <Td mono className="sticky left-0 bg-gray-50 z-10">
                column sum (p&lt;α ∧ z&gt;0)
              </Td>
              {colSupportSums.map((s, j) => (
                <Td key={`sum-${j}`} center className="font-mono">{s}</Td>
              ))}
            </tr>
          </tfoot>
        </table>
      </div>
    </section>
  );
}

/** ---------- Agreement Count Matrix (raw + normalized column sums) + CSV ---------- */
function CountsMatrixTable({
  title,
  nodes,
  sizes,
  counts,
  maxCount,
  jobId,
  method,
  metric,
  alpha,
}: {
  title: string;
  nodes: string[];
  sizes: Record<string, number>;
  counts: number[][];
  maxCount: number; // number of available methods (dynamic)
  jobId: string;
  method: MethodName;
  metric: "final" | "init";
  alpha: number;
}) {
  const includeDiag = false;
  const N = nodes.length;

  const colSums = Array.from({ length: N }, (_, j) =>
    counts.reduce((acc, row, i) => {
      const val = row?.[j];
      if (!Number.isFinite(val)) return acc;
      if (!includeDiag && i === j) return acc;
      return acc + (val as number);
    }, 0)
  );

  // Normalize by (N-1)*maxCount; label dynamically
  const denom = (includeDiag ? N : N - 1) * maxCount;
  const denomLabel = `${includeDiag ? N : N - 1}×${maxCount}`;
  const colNorms = colSums.map((s) => (Number.isFinite(s) && denom > 0 ? s / denom : NaN));

  const fmtSig2 = (x: number) =>
    Number.isFinite(x)
      ? new Intl.NumberFormat(undefined, { minimumSignificantDigits: 2, maximumSignificantDigits: 2 }).format(x)
      : "—";

  function cellClass(c: number, i: number, j: number): string {
    if (i === j) return "bg-gray-100 text-gray-400";
    if (!Number.isFinite(c)) return "";
    if (c >= maxCount) return "bg-emerald-300";
    if (c === Math.max(1, maxCount - 1)) return "bg-emerald-200";
    if (c >= 1) return "bg-emerald-100";
    return "bg-white";
  }

  const downloadCountsCsv = () => {
    const rows: (string | number)[][] = [
      ["row/col", ...nodes.map(showNode)],
      ["sample size (per node)", ...nodes.map((n) => sizes[n] ?? 0)],
    ];
    for (let i = 0; i < N; i++) {
      const r: (string | number)[] = [showNode(nodes[i])];
      for (let j = 0; j < N; j++) {
        r.push(i === j ? "—" : (counts[i]?.[j] ?? ""));
      }
      rows.push(r);
    }
    rows.push(["column sum (excl. diag)", ...colSums]);
    rows.push([`column sum ÷ (${denomLabel})`, ...colNorms.map(v => (Number.isFinite(v) ? Number(v.toPrecision(2)) : "—"))]);

    const fname = `wmw_agree_${metric}_${method}_${jobId}_alpha-${alpha}.csv`;
    downloadCsv(fname, rows);
  };

  return (
    <section className="space-y-2">
      <div className="flex items-center justify-between">
        <h2 className="text-lg font-semibold">{title}</h2>
        <button
          type="button"
          onClick={downloadCountsCsv}
          className="px-2 py-1 text-xs rounded bg-black text-white"
          title="Download counts + footers CSV"
        >
          CSV (counts + footers)
        </button>
      </div>

      <div className="overflow-auto border rounded-xl">
        <table className="min-w-max text-sm">
          <thead className="sticky top-0 bg-white z-10">
            <tr>
              <Th className="sticky left-0 bg-white z-20">row ↓ / col →</Th>
              {nodes.map((n) => (
                <Th key={n} className="text-center">{showNode(n)}</Th>
              ))}
            </tr>
            <tr>
              <Th className="sticky left-0 bg-white z-20">sample size (per node)</Th>
              {nodes.map((n) => (
                <Th key={`${n}-n`} className="text-center">{sizes[n] ?? 0}</Th>
              ))}
            </tr>
          </thead>

          <tbody>
            {nodes.map((rowName, i) => (
              <tr key={rowName} className="hover:bg-gray-50">
                <Td mono className="sticky left-0 bg-white z-10">{showNode(rowName)}</Td>
                {nodes.map((_, j) => (
                  <Td key={`${i}-${j}`} center className={cellClass(counts[i]?.[j], i, j)}>
                    {i === j ? "—" : <span className="font-mono">{counts[i]?.[j] ?? ""}</span>}
                  </Td>
                ))}
              </tr>
            ))}
          </tbody>

          <tfoot className="bg-gray-50">
            <tr className="font-semibold">
              <Td mono className="sticky left-0 bg-gray-50 z-10">
                column sum{includeDiag ? "" : " (excl. diag)"}
              </Td>
              {colSums.map((s, j) => (
                <Td key={`sum-${j}`} center className="font-mono">{s}</Td>
              ))}
            </tr>
            <tr className="font-semibold">
              <Td mono className="sticky left-0 bg-gray-50 z-10">
                column sum ÷ ({denomLabel})
              </Td>
              {colNorms.map((v, j) => (
                <Td key={`norm-${j}`} center className="font-mono">
                  {fmtSig2(v)}
                </Td>
              ))}
            </tr>
          </tfoot>
        </table>
      </div>
      <p className="text-xs text-gray-500">
        Raw sums add agreement counts down each column (excluding the diagonal).
        Normalized values divide by {denomLabel} (max attainable per column with the available methods).
      </p>
    </section>
  );
}

/** ---------- Legends & UI atoms ---------- */

function Legend({ alpha }: { alpha: number }) {
  return (
    <div className="flex items-center gap-4 text-xs text-gray-600">
      <div className="flex items-center gap-2">
        <span className="inline-block w-4 h-4 rounded bg-emerald-100 border" />
        <span>p &lt; α and z &gt; 0 (evidence column &gt; row)</span>
      </div>
      <div className="flex items-center gap-2">
        <span className="inline-block w-4 h-4 rounded bg-rose-50 border" />
        <span>medians favor row ≥ column</span>
      </div>
      <div>α = {alpha}</div>
    </div>
  );
}

function AgreementLegend() {
  return (
    <div className="flex items-center gap-4 text-xs text-gray-600">
      <div className="flex items-center gap-2">
        <span className="inline-block w-4 h-4 rounded bg-emerald-100 border" />
        <span>1 method supports col &gt; row</span>
      </div>
      <div className="flex items-center gap-2">
        <span className="inline-block w-4 h-4 rounded bg-emerald-200 border" />
        <span>2 methods support col &gt; row</span>
      </div>
      <div className="flex items-center gap-2">
        <span className="inline-block w-4 h-4 rounded bg-emerald-300 border" />
        <span>All available methods support col &gt; row</span>
      </div>
    </div>
  );
}

function Th({
  children,
  className,
  colSpan,
  rowSpan,
}: {
  children: React.ReactNode;
  className?: string;
  colSpan?: number;
  rowSpan?: number;
}) {
  return (
    <th
      className={`text-left font-semibold px-3 py-2 whitespace-nowrap border-b ${className ?? ""}`}
      colSpan={colSpan}
      rowSpan={rowSpan}
    >
      {children}
    </th>
  );
}

function Td({
  children,
  mono,
  num,
  center,
  className,
}: {
  children: React.ReactNode;
  mono?: boolean;
  num?: boolean;
  center?: boolean;
  className?: string;
}) {
  const align = center ? "text-center" : num ? "text-right" : "text-left";
  const numCls = num ? "tabular-nums" : "";
  const monoCls = mono ? "font-mono" : "";
  return (
    <td className={`px-3 py-2 ${align} ${numCls} ${monoCls} border-b ${className ?? ""}`}>
      {children}
    </td>
  );
}
