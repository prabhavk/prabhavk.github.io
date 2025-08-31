// app/wmw_comp_method/page.tsx
"use client";

const showNode = (s: string) => s.replace(/^h_/, "");

import React, { useEffect, useMemo, useState } from "react";

/* ------------ Types (aligned with /api/wmw_comp_method) ------------ */

type MethodName = "Parsimony" | "Dirichlet" | "SSH";
type PairKey = "parsimony_vs_dirichlet" | "parsimony_vs_ssh" | "dirichlet_vs_ssh";

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
  pairs: Record<PairKey, PairResult>;
};

/* ----------------------- Small Utils ----------------------- */

const NODES: string[] = Array.from({ length: 17 }, (_, i) => `h_${21 + i}`);
const INIT_OPTIONS = ["(all nodes)", ...NODES] as const;

function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isApiErr(x: unknown): x is { error: string } {
  return isRecord(x) && typeof x["error"] === "string";
}

async function fetchJSON<T>(url: string, init?: RequestInit): Promise<T> {
  const res = await fetch(url, { ...init, headers: { accept: "application/json", ...(init?.headers || {}) } });
  const ct = res.headers.get("content-type") || "";
  const body = await res.text(); // read once for better diagnostics

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

/* ----------------------- Page Component ----------------------- */

type TableRow = {
  node: string;
  sizes: Record<MethodName, number>;
  // Parsimony vs Dirichlet
  pdZ: number; pdP: number; pdWinner: MethodName | "none";
  // Parsimony vs SSH
  psZ: number; psP: number; psWinner: MethodName | "none";
  // Dirichlet vs SSH
  dsZ: number; dsP: number; dsWinner: MethodName | "none";
};

export default function WmwCompMethodPage() {
  const [jobId, setJobId] = useState<string>("");
  const [alpha, setAlpha] = useState<number>(0.05);
  const [initChoice, setInitChoice] = useState<(typeof INIT_OPTIONS)[number]>("(all nodes)");

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

  // Auto-compute:
  // - on first load (once jobId is available),
  // - whenever initChoice (initialization criterion) changes,
  // - whenever alpha changes.
  useEffect(() => {
    if (!jobId) return;
    void compute();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [jobId, initChoice, alpha]);

  async function compute() {
    if (!jobId) {
      setErr("No job selected. Pick a job on the Precomputed Results page.");
      return;
    }
    setLoading(true);
    setErr(null);
    try {
      const selectedNodes = initChoice === "(all nodes)" ? NODES : [initChoice];
      const results = await Promise.all(
        selectedNodes.map(async (node) => {
          const u = new URL("/api/wmw_comp_method", window.location.origin);
          u.searchParams.set("job", jobId);
          u.searchParams.set("alpha", String(alpha));
          u.searchParams.set("node", node);

          const json = await fetchJSON<ApiResp>(u.toString(), { cache: "no-store" });

          const pd = json.pairs.parsimony_vs_dirichlet;
          const ps = json.pairs.parsimony_vs_ssh;
          const ds = json.pairs.dirichlet_vs_ssh;

          const row: TableRow = {
            node,
            sizes: json.sizes,
            pdZ: pd.z, pdP: pd.p_value, pdWinner: pd.winner,
            psZ: ps.z, psP: ps.p_value, psWinner: ps.winner,
            dsZ: ds.z, dsP: ds.p_value, dsWinner: ds.winner,
          };
          return row;
        })
      );
      setRows(results);
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : "Failed to compute";
      setErr(msg);
      setRows([]);
    } finally {
      setLoading(false);
    }
  }

  // CSV download of the visible table
  function downloadCSV() {
  const header = [
    "node",
    "n_Parsimony",
    "n_Dirichlet",
    "n_SSH",
    "n_total",
    "pd_z", "pd_p", "pd_winner",
    "ps_z", "ps_p", "ps_winner",
    "ds_z", "ds_p", "ds_winner",
  ];

  const lines = rows.map((r) => {
    const nPar = r.sizes.Parsimony ?? 0;
    const nDir = r.sizes.Dirichlet ?? 0;
    const nSSH = r.sizes.SSH ?? 0;
    const nTot = nPar + nDir + nSSH;

    return [
      showNode(showNode(r.node)),
      String(nPar),
      String(nDir),
      String(nSSH),
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
    a.download = `wmw_comp_method_all-nodes_alpha-${alpha}.csv`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }

  function safeCSV(cell: string): string {
    // Wrap if contains commas/quotes/newlines
    if (/[",\n]/.test(cell)) {
      return `"${cell.replace(/"/g, '""')}"`;
    }
    return cell;
  }

  const alphaStr = useMemo(() => alpha.toString(), [alpha]);

  return (
    <div className="p-6 max-w-[95vw] mx-auto space-y-6">
      <header className="space-y-1">
        <h1 className="text-2xl font-bold">Method Comparison (WMW one-sided)</h1>
        <p className="text-sm text-gray-600">
          For each selected node, we compare methods pairwise. Cells show one-tailed p for H₁: method with larger median wins.
        </p>
      </header>

      <section className="flex flex-wrap items-end gap-3">
        <label className="flex flex-col">
          <span className="text-sm">Initialization criterion (node/root)</span>
          <select
            value={initChoice}
            onChange={(e) => setInitChoice(e.target.value as (typeof INIT_OPTIONS)[number])}
            className="border rounded-md px-3 py-2"
          >
            {INIT_OPTIONS.map((opt) => <option key={opt} value={opt}>{opt}</option>)}
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
          />
        </label>

        <div className="text-sm text-gray-600">
          Job: <span className="font-mono">{jobId || "(none)"} </span>
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
          <thead className="bg-white">
            <tr>
              <Th rowSpan={2}>root</Th>
              <Th colSpan={4} className="text-center">sample size</Th>
              <Th colSpan={3} className="text-center">Parsimony vs Dirichlet</Th>
              <Th colSpan={3} className="text-center">Parsimony vs SSH</Th>
              <Th colSpan={3} className="text-center">Dirichlet vs SSH</Th>
            </tr>
            <tr>
              <Th className="text-center">Parsimony</Th>
              <Th className="text-center">Dirichlet</Th>
              <Th className="text-center">SSH</Th>
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
                const nSSH = r.sizes.SSH ?? 0;
                return (
                  <tr key={showNode(r.node)} className="hover:bg-gray-50">
                    <Td mono>{showNode(r.node)}</Td>
                    <Td num center>{nPar}</Td>
                    <Td num center>{nDir}</Td>
                    <Td num center>{nSSH}</Td>
                    <Td num center>{nPar + nDir + nSSH}</Td>

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
        Table auto-updates when α or the initialization criterion changes. “winner” is the method with larger median when the one-sided test is significant.
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
