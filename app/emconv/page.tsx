// app/iters/page.tsx
"use client";

import React, { useEffect, useMemo, useState } from "react";

type MethodName = "Parsimony" | "Dirichlet" | "SSH";

type Stats = { n: number; avg: number; max: number };
type ApiResp = {
  job_id: string;
  roots: string[];
  stats: Record<string, Record<MethodName, Stats>>;
};
type ApiErr = { error: string };

type TableRow = {
  root: string;
  P_avg?: number; P_max?: number;
  D_avg?: number; D_max?: number;
  S_avg?: number; S_max?: number;
};

function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isErr(x: unknown): x is ApiErr {
  return isRecord(x) && typeof x["error"] === "string";
}

function fmtAvg(v?: number) {
  return Number.isFinite(v!) ? (v as number).toFixed(1) : "";
}
function fmtMax(v?: number) {
  return Number.isFinite(v!) ? String(Math.round(v as number)) : "";
}

export default function ItersPage() {
  const [job, setJob] = useState<string>("");
  const [rows, setRows] = useState<TableRow[]>([]);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  useEffect(() => {
    try {
      const saved = localStorage.getItem("emtr:selectedJobId") ?? "";
      setJob(saved);
    } catch { /* ignore */ }
  }, []);

  async function load() {
    if (!job) {
      setErr("No job selected. Please pick a job on the Precomputed Results page.");
      setRows([]);
      return;
    }
    setLoading(true);
    setErr(null);
    try {
      const u = new URL("/api/iters", window.location.origin);
      u.searchParams.set("job", job);
      const res = await fetch(u.toString(), { cache: "no-store" });

      const ct = res.headers.get("content-type") || "";
      if (!ct.includes("application/json")) {
        const body = await res.text().catch(() => "");
        throw new Error(`Expected JSON, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0,160)}`);
      }

      const json: unknown = await res.json();
      if (!res.ok || isErr(json)) throw new Error(isErr(json) ? json.error : `HTTP ${res.status}`);
      const data = json as ApiResp;

      const tbl: TableRow[] = data.roots.map((root) => {
        const m = data.stats[root] || {};
        const P = m["Parsimony"]; const D = m["Dirichlet"]; const S = m["SSH"];
        return {
          root,
          P_avg: P?.avg, P_max: P?.max,
          D_avg: D?.avg, D_max: D?.max,
          S_avg: S?.avg, S_max: S?.max,
        };
      });
      setRows(tbl);
    } catch (e: unknown) {
      setErr(e instanceof Error ? e.message : "Failed to load");
      setRows([]);
    } finally {
      setLoading(false);
    }
  }

  const hasData = useMemo(() => rows.length > 0, [rows]);

  return (
    <div className="p-6 max-w-6xl mx-auto space-y-6">
      <header className="flex items-end gap-3">
        <h1 className="text-2xl font-bold">EM iterations</h1>
        <div className="text-sm text-gray-600">
          {job ? <>Job: <span className="font-mono">{job}</span></> : "No job selected"}
        </div>
        <button
          className="ml-auto px-4 py-2 rounded bg-black text-white disabled:opacity-50"
          onClick={() => void load()}
          disabled={loading}
        >
          {loading ? "Loading…" : "Compute"}
        </button>
      </header>

      <div className="overflow-auto border rounded-xl">
        <table className="min-w-full text-sm">
          <thead className="bg-black text-white">
            <tr>
              <Th rowSpan={2}>root</Th>
              <Th colSpan={2} className="text-center">Parsimony</Th>
              <Th colSpan={2} className="text-center">Dirichlet</Th>
              <Th colSpan={2} className="text-center">SSH</Th>
            </tr>
            <tr>
              <Th className="text-center">avg iters</Th>
              <Th className="text-center">max iters</Th>
              <Th className="text-center">avg iters</Th>
              <Th className="text-center">max iters</Th>
              <Th className="text-center">avg iters</Th>
              <Th className="text-center">max iters</Th>
            </tr>
          </thead>
          <tbody>
            {err ? (
              <tr><td colSpan={7} className="p-6 text-center text-red-600">{err}</td></tr>
            ) : !hasData ? (
              <tr>
                <td colSpan={7} className="p-6 text-center">
                  No results. {job ? "Click Compute." : "Select a job on the Precomputed Results page."}
                </td>
              </tr>
            ) : (
              rows.map((r) => (
                <tr key={r.root} className="hover:bg-gray-50">
                  <Td mono>{r.root}</Td>
                  <Td center num>{fmtAvg(r.P_avg)}</Td>
                  <Td center num>{fmtMax(r.P_max)}</Td>
                  <Td center num>{fmtAvg(r.D_avg)}</Td>
                  <Td center num>{fmtMax(r.D_max)}</Td>
                  <Td center num>{fmtAvg(r.S_avg)}</Td>
                  <Td center num>{fmtMax(r.S_max)}</Td>
                </tr>
              ))
            )}
          </tbody>
        </table>
      </div>

      <p className="text-xs text-gray-500">
        Averages and maxima are computed across repetitions for each (root, method).
      </p>
    </div>
  );
}

function Th({
  children, className, colSpan, rowSpan,
}: { children: React.ReactNode; className?: string; colSpan?: number; rowSpan?: number }) {
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
  children, mono, num, center, className,
}: { children: React.ReactNode; mono?: boolean; num?: boolean; center?: boolean; className?: string }) {
  const align = center ? "text-center" : num ? "text-right" : "text-left";
  const numCls = num ? "tabular-nums" : "";
  const monoCls = mono ? "font-mono" : "";
  return (
    <td className={`px-3 py-2 ${align} ${numCls} ${monoCls} ${className ?? ""}`}>
      {children}
    </td>
  );
}
