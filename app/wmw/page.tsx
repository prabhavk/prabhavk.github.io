// app/wmw/page.tsx
"use client";

import React, { useEffect, useMemo, useState } from "react";

type MethodName = "Parsimony" | "Dirichlet" | "SSH";
type PairKey = "parsimony_vs_dirichlet" | "parsimony_vs_ssh" | "dirichlet_vs_ssh";

type PairResult = {
  winner: MethodName | "none";
  p_value: number;
  z: number;
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

type TableRow = {
  node: string;
  // Parsimony vs Dirichlet
  pdZ: number;
  pdP: number;
  pdWinner: MethodName | "none";
  // Parsimony vs SSH
  psZ: number;
  psP: number;
  psWinner: MethodName | "none";
  // Dirichlet vs SSH
  dsZ: number;
  dsP: number;
  dsWinner: MethodName | "none";
};

const NODES: string[] = Array.from({ length: 17 }, (_, i) => `h_${21 + i}`);

export default function WmwPage() {
  const [jobId, setJobId] = useState<string>("");
  const [alpha, setAlpha] = useState<number>(0.05);
  const [rows, setRows] = useState<TableRow[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId");
      if (saved) setJobId(saved);
    } catch {
      /* ignore */
    }
  }, []);

  async function computeAll() {
    if (!jobId) {
      setError("No job selected. Pick a job on the Precomputed Results page.");
      return;
    }
    setLoading(true);
    setError(null);
    try {
      const results = await Promise.all(
        NODES.map(async (node) => {
          const u = new URL("/api/mwu", window.location.origin);
          u.searchParams.set("job", jobId);
          u.searchParams.set("node", node);
          u.searchParams.set("alpha", String(alpha));
          const res = await fetch(u.toString(), { cache: "no-store" });
          const json: unknown = await res.json();
          if (!res.ok) {
            const msg = isApiErr(json) ? json.error : `HTTP ${res.status}`;
            throw new Error(msg);
          }
          if (!isApiResp(json)) throw new Error("Invalid response from /api/mwu");

          const pd = json.pairs.parsimony_vs_dirichlet;
          const ps = json.pairs.parsimony_vs_ssh;
          const ds = json.pairs.dirichlet_vs_ssh;

          const row: TableRow = {
            node,
            pdZ: pd.z,
            pdP: pd.p_value,
            pdWinner: pd.winner,
            psZ: ps.z,
            psP: ps.p_value,
            psWinner: ps.winner,
            dsZ: ds.z,
            dsP: ds.p_value,
            dsWinner: ds.winner,
          };
          return row;
        })
      );
      setRows(results);
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : "Failed to compute";
      setError(msg);
      setRows([]);
    } finally {
      setLoading(false);
    }
  }

  const alphaStr = useMemo(() => alpha.toString(), [alpha]);

  return (
    <div className="p-6 max-w-6xl mx-auto">
      <h1 className="text-2xl font-bold mb-3">Wilcoxon–Mann–Whitney (one-sided)</h1>
      <div className="text-sm text-gray-600 mb-4">
        Selected job: <span className="font-mono">{jobId || "(none)"}</span>
      </div>

      <div className="flex items-end gap-3 mb-4">
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

        <button
          type="button"
          onClick={() => void computeAll()}
          className="px-4 py-2 rounded-md bg-black text-white disabled:opacity-50"
          disabled={!jobId || loading}
        >
          {loading ? "Computing…" : "Compute"}
        </button>
      </div>

      <div className="overflow-auto border rounded-xl">
        <table className="min-w-full text-sm">
          <thead className="bg-black text-white">
            <tr>
              <Th rowSpan={2}>root</Th>
              <Th colSpan={3} className="text-center">Parsimony vs Dirichlet</Th>
              <Th colSpan={3} className="text-center">Parsimony vs SSH</Th>
              <Th colSpan={3} className="text-center">Dirichlet vs SSH</Th>
            </tr>
            <tr>
              {/* PD */}
              <Th className="text-center">z-score</Th>
              <Th className="text-center">p-value</Th>
              <Th className="text-center">winner</Th>
              {/* PS */}
              <Th className="text-center">z-score</Th>
              <Th className="text-center">p-value</Th>
              <Th className="text-center">winner</Th>
              {/* DS */}
              <Th className="text-center">z-score</Th>
              <Th className="text-center">p-value</Th>
              <Th className="text-center">winner</Th>
            </tr>
          </thead>
          <tbody>
            {error ? (
              <tr>
                <td colSpan={10} className="p-6 text-center text-red-600">{error}</td>
              </tr>
            ) : rows.length === 0 ? (
              <tr>
                <td colSpan={10} className="p-6 text-center">No results yet. Click <em>Compute</em>.</td>
              </tr>
            ) : (
              rows.map((r) => (
                <tr key={r.node} className="hover:bg-gray-50">
                  <Td mono>{r.node}</Td>

                  {/* PD group: z, p, winner */}
                  <Td num center>{fmtZ(r.pdZ)}</Td>
                  <Td num center>{fmtP(r.pdP)}</Td>
                  <Td center>{fmtWinner(r.pdWinner)}</Td>

                  {/* PS group: z, p, winner */}
                  <Td num center>{fmtZ(r.psZ)}</Td>
                  <Td num center>{fmtP(r.psP)}</Td>
                  <Td center>{fmtWinner(r.psWinner)}</Td>

                  {/* DS group: z, p, winner */}
                  <Td num center>{fmtZ(r.dsZ)}</Td>
                  <Td num center>{fmtP(r.dsP)}</Td>
                  <Td center>{fmtWinner(r.dsWinner)}</Td>
                </tr>
              ))
            )}
          </tbody>
        </table>
      </div>

      <p className="text-xs text-gray-500 mt-4">
        Winner is the method with larger median log-likelihood when the one-sided test is significant (p = 1 − Φ(|z|)).
      </p>
    </div>
  );
}

// ---------- helpers & type guards ----------

function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isPairResult(x: unknown): x is PairResult {
  return (
    isRecord(x) &&
    (x.winner === "Parsimony" || x.winner === "Dirichlet" || x.winner === "SSH" || x.winner === "none") &&
    typeof x.p_value === "number" &&
    typeof x.z === "number" &&
    typeof x.nA === "number" &&
    typeof x.nB === "number"
  );
}
function isApiErr(x: unknown): x is { error: string } {
  return isRecord(x) && typeof x["error"] === "string";
}
function isApiResp(x: unknown): x is ApiResp {
  if (!isRecord(x) || !isRecord(x["pairs"])) return false;
  const p = x["pairs"] as Record<string, unknown>;
  return (
    isPairResult(p["parsimony_vs_dirichlet"]) &&
    isPairResult(p["parsimony_vs_ssh"]) &&
    isPairResult(p["dirichlet_vs_ssh"])
  );
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
      className={`text-left font-semibold px-3 py-2 whitespace-nowrap ${className ?? ""}`}
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
    <td className={`px-3 py-2 ${align} ${numCls} ${monoCls} ${className ?? ""}`}>
      {children}
    </td>
  );
}
