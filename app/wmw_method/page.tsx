// app/wmw_comp_method/page.tsx
"use client";

const showNode = (s: string) => s.replace(/^h_/, "");

import React, { useEffect, useMemo, useState } from "react";

/* ------------ Types (aligned with /api/wmw_comp_method) ------------ */

type MethodName = "Parsimony" | "Dirichlet" | "HSS";
type PairKey = "parsimony_vs_dirichlet" | "parsimony_vs_hss" | "dirichlet_vs_hss";

type PairResult = {
  winner: MethodName | "none";
  p_value: number;   // one-sided p (winner's direction)
  z: number;         // z-statistic from MWU
  nA: number;
  nB: number;
};

type ApiResp = {
  job_id: string;
  node?: string;
  alpha: number;
  sizes: Record<MethodName, number>;
  pairs: Record<PairKey, PairResult> | Record<string, PairResult>; // tolerate legacy keys
};

type RootsResp =
  | { ok: true; roots: string[] }
  | { ok: false; error: string };

/* ----------------------- Small Utils ----------------------- */

function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isApiErr(x: unknown): x is { error: string } {
  return isRecord(x) && typeof x["error"] === "string";
}

async function fetchJSON<T>(url: string, init?: RequestInit): Promise<T> {
  const res = await fetch(url, { ...init, headers: { accept: "application/json", ...(init?.headers || {}) } });
  const ct = res.headers.get("content-type") || "";
  const body = await res.text();

  if (!ct.toLowerCase().includes("application/json")) {
    throw new Error(`Expected JSON but got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 200)}`);
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

function fmtP(p: number): string {
  if (!Number.isFinite(p)) return "";
  if (p === 0) return "0";
  if (p < 1e-3) return p.toExponential(2);
  return p.toFixed(4);
}
function fmtZ(z: number): string {
  if (!Number.isFinite(z)) return "";
  return z.toFixed(3);
}
function fmtWinner(w: MethodName | "none"): string {
  return w === "none" ? "" : w;
}

function numericNodeId(s: string): number {
  const m = /^h_(\d+)$/.exec(s);
  return m ? Number(m[1]) : Number.POSITIVE_INFINITY;
}

/* ---------- Type guards for pairs (to avoid `any`) ---------- */

function isMethodNameOrNone(x: unknown): x is MethodName | "none" {
  return x === "Parsimony" || x === "Dirichlet" || x === "HSS" || x === "none";
}

function isPairResult(x: unknown): x is PairResult {
  return (
    isRecord(x) &&
    isMethodNameOrNone(x.winner) &&
    typeof x.p_value === "number" &&
    typeof x.z === "number" &&
    typeof x.nA === "number" &&
    typeof x.nB === "number"
  );
}

/** Safely grab a pair by possible keys (supports legacy aliases). */
function pickPair(pairs: ApiResp["pairs"], keys: string[]): PairResult | undefined {
  if (!isRecord(pairs)) return undefined;
  for (const k of keys) {
    const v = (pairs as Record<string, unknown>)[k];
    if (isPairResult(v)) return v;
  }
  return undefined;
}

/* ----------------------- Page Component ----------------------- */

type TableRow = {
  node: string;
  sizes: Record<MethodName, number>;
  // Parsimony vs Dirichlet
  pdZ: number; pdP: number; pdWinner: MethodName | "none";
  // Parsimony vs HSS
  psZ: number; psP: number; psWinner: MethodName | "none";
  // Dirichlet vs HSS
  dsZ: number; dsP: number; dsWinner: MethodName | "none";
};

export default function WmwCompMethodPage() {
  const [jobId, setJobId] = useState<string>("");
  const [alpha, setAlpha] = useState<number>(0.05);

  // Dynamic roots from the server (DB)
  const [rootOptions, setRootOptions] = useState<string[]>([]);

  const [rows, setRows] = useState<TableRow[]>([]);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  // Load job id from localStorage, once.
  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId");
      if (saved) setJobId(saved);
    } catch {/* ignore */}
  }, []);

  // Load available roots whenever jobId changes
  useEffect(() => {
    if (!jobId) return;

    (async () => {
      try {
        const u = new URL("/api/wmw_comp_method/roots", window.location.origin);
        u.searchParams.set("job", jobId);
        const j = await fetchJSON<RootsResp>(u.toString(), { cache: "no-store" });

        if (j.ok) {
          const uniqueSorted = Array.from(new Set(j.roots || []))
            .sort((a, b) => numericNodeId(a) - numericNodeId(b));
          setRootOptions(uniqueSorted);
        } else {
          setRootOptions([]);
        }
      } catch {
        setRootOptions([]);
      }
    })();
  }, [jobId]);

  // Auto-compute when job, root list, or alpha change
  useEffect(() => {
    if (!jobId) return;
    if (!rootOptions.length) {
      setRows([]);
      return;
    }
    void computeAllRoots();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [jobId, rootOptions, alpha]);

  async function computeAllRoots() {
    if (!jobId) {
      setErr("No job selected. Pick a job on the Precomputed Results page.");
      return;
    }
    if (!rootOptions.length) {
      setErr("No roots found for this job.");
      setRows([]);
      return;
    }

    setLoading(true);
    setErr(null);
    try {
      const results = await Promise.all(
        rootOptions.map(async (node) => {
          const u = new URL("/api/wmw_comp_method", window.location.origin);
          u.searchParams.set("job", jobId);
          u.searchParams.set("alpha", String(alpha));
          u.searchParams.set("node", node); // exact DB value

          const json = await fetchJSON<ApiResp>(u.toString(), { cache: "no-store" });

          const pd = pickPair(json.pairs, ["parsimony_vs_dirichlet"]);
          // tolerate either ..._hss or older ..._ssh keys from backend
          const ps = pickPair(json.pairs, ["parsimony_vs_hss", "parsimony_vs_ssh"]);
          const ds = pickPair(json.pairs, ["dirichlet_vs_hss", "dirichlet_vs_ssh"]);

          const row: TableRow = {
            node,
            sizes: json.sizes,
            pdZ: pd?.z ?? NaN, pdP: pd?.p_value ?? NaN, pdWinner: pd?.winner ?? "none",
            psZ: ps?.z ?? NaN, psP: ps?.p_value ?? NaN, psWinner: ps?.winner ?? "none",
            dsZ: ds?.z ?? NaN, dsP: ds?.p_value ?? NaN, dsWinner: ds?.winner ?? "none",
          };
          return row;
        })
      );
      results.sort((a, b) => numericNodeId(a.node) - numericNodeId(b.node));
      setRows(results);
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : "Failed to compute";
      setErr(msg);
      setRows([]);
    } finally {
      setLoading(false);
    }
  }

  // CSV download of the visible table (all roots)
  function downloadCSV() {
    const header = [
      "root",
      "n_Parsimony",
      "n_Dirichlet",
      "n_HSS",
      "n_total",
      "pd_z", "pd_p", "pd_winner",
      "ps_z", "ps_p", "ps_winner",
      "ds_z", "ds_p", "ds_winner",
    ];

    const lines = rows.map((r) => {
      const nPar = r.sizes.Parsimony ?? 0;
      const nDir = r.sizes.Dirichlet ?? 0;
      const nHSS = r.sizes.HSS ?? 0;
      const nTot = nPar + nDir + nHSS;

      return [
        showNode(r.node),
        String(nPar),
        String(nDir),
        String(nHSS),
        String(nTot),
        fmtZ(r.pdZ), fmtP(r.pdP), fmtWinner(r.pdWinner),
        fmtZ(r.psZ), fmtP(r.psP), fmtWinner(r.psWinner),
        fmtZ(r.dsZ), fmtP(r.dsP), fmtWinner(r.dsWinner),
      ];
    });

    const csv = [header.join(","), ...lines.map(cols => cols.map(safeCSV).join(","))].join("\n");
    const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `wmw_comp_method_all-roots_alpha-${alpha}.csv`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }

  function safeCSV(cell: string): string {
    if (/[",\n]/.test(cell)) return `"${cell.replace(/"/g, '""')}"`;
    return cell;
  }

  const alphaStr = useMemo(() => alpha.toString(), [alpha]);

  return (
    <div className="p-6 max-w-[95vw] mx-auto space-y-6">
      <header className="space-y-1">
        <h1 className="text-2xl text-black font-bold">Method Comparison (WMW one-sided)</h1>
        <p className="text-sm text-black">
          One row per root. Cells show one-tailed p for H₁: method with larger median wins.
        </p>
      </header>

      <section className="flex flex-wrap items-end gap-3">
        <label className="flex flex-col">
          <span className="text-black text-sm">α (significance)</span>
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
          />
        </label>

        <div className="text-sm text-gray-600">
          Job: <span className="font-mono">{jobId || "(none)"} </span>
          {rootOptions.length ? (
            <span className="ml-2 text-gray-500">roots: {rootOptions.length}</span>
          ) : null}
        </div>

        <div className="ml-auto flex items-center gap-3">
          <button
            type="button"
            onClick={downloadCSV}
            className="px-3 py-2 rounded-md bg-black text-white disabled:opacity-50 text-sm"
            disabled={rows.length === 0 || !!err}
            title="Download the visible table as CSV"
          >
            Download CSV
          </button>
        </div>
      </section>

      {err && (
        <div className="p-3 rounded-md bg-red-50 text-red-700 border border-red-200">
          {err}
        </div>
      )}

      <div className="overflow-auto border rounded-xl">
        <table className="min-w-full text-sm">
          <thead className="bg-white text-black">
            <tr>
              <Th rowSpan={2}>root</Th>
              <Th colSpan={4} className="text-center">sample size</Th>
              <Th colSpan={3} className="text-center">Parsimony vs Dirichlet</Th>
              <Th colSpan={3} className="text-center">Parsimony vs HSS</Th>
              <Th colSpan={3} className="text-center">Dirichlet vs HSS</Th>
            </tr>
            <tr>
              <Th className="text-center">Parsimony</Th>
              <Th className="text-center">Dirichlet</Th>
              <Th className="text-center">HSS</Th>
              <Th className="text-center">total</Th>

              <Th className="text-center">z</Th>
              <Th className="text-center">p</Th>
              <Th className="text-center">winner</Th>

              <Th className="text-center">z</Th>
              <Th className="text-center">p</Th>
              <Th className="text-center">winner</Th>

              <Th className="text-center">z</Th>
              <Th className="text-center">p</Th>
              <Th className="text-center">winner</Th>
            </tr>
          </thead>
          <tbody>
            {loading && rows.length === 0 ? (
              <tr><Td colSpan={14} center>Computing…</Td></tr>
            ) : rows.length === 0 ? (
              <tr><Td colSpan={14} center>No results to display.</Td></tr>
            ) : (
              rows.map((r) => {
                const nPar = r.sizes.Parsimony ?? 0;
                const nDir = r.sizes.Dirichlet ?? 0;
                const nHSS = r.sizes.HSS ?? 0;
                return (
                  <tr key={r.node} className="text-black hover:bg-gray-50">
                    <Td mono>{showNode(r.node)}</Td>
                    <Td num center>{nPar}</Td>
                    <Td num center>{nDir}</Td>
                    <Td num center>{nHSS}</Td>
                    <Td num center>{nPar + nDir + nHSS}</Td>

                    <Td num center>{fmtZ(r.pdZ)}</Td>
                    <Td num center>{fmtP(r.pdP)}</Td>
                    <Td center>{fmtWinner(r.pdWinner)}</Td>

                    <Td num center>{fmtZ(r.psZ)}</Td>
                    <Td num center>{fmtP(r.psP)}</Td>
                    <Td center>{fmtWinner(r.psWinner)}</Td>

                    <Td num center>{fmtZ(r.dsZ)}</Td>
                    <Td num center>{fmtP(r.dsP)}</Td>
                    <Td center>{fmtWinner(r.dsWinner)}</Td>
                  </tr>
                );
              })
            )}
          </tbody>
        </table>
      </div>

      <p className="text-xs text-gray-500">
        The table lists one row per root. It refreshes when α changes.
      </p>
    </div>
  );
}

/* ----------------------- Table atoms ----------------------- */

function Th({
  children, className, colSpan, rowSpan,
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
  children, mono, num, center, className, colSpan,
}: {
  children: React.ReactNode;
  mono?: boolean;
  num?: boolean;
  center?: boolean;
  className?: string;
  colSpan?: number;
}) {
  const align = center ? "text-center" : num ? "text-right" : "text-left";
  const numCls = num ? "tabular-nums" : "";
  const monoCls = mono ? "font-mono" : "";
  return (
    <td
      className={`px-3 py-2 ${align} ${numCls} ${monoCls} border-b ${className ?? ""}`}
      colSpan={colSpan}
    >
      {children}
    </td>
  );
}
