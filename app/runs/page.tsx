// app/runs/page.tsx
"use client";
import React, { useEffect, useMemo, useState } from "react";

type RunRow = {
  job_id: string;
  status: string;
  created_at: string;
  finished_at: string;
  dur_minutes: number | null;
  thr: number | null;
  reps: number | null;
  max_iter: number | null;

  d_pi_1: number | null;
  d_pi_2: number | null;
  d_pi_3: number | null;
  d_pi_4: number | null;
  d_m_1: number | null;
  d_m_2: number | null;
  d_m_3: number | null;
  d_m_4: number | null;
};

type ApiOk = { page: number; pageSize: number; total: number; rows: RunRow[] };
type ApiErr = { error: string };

function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}
function isApiErr(x: unknown): x is ApiErr {
  return isRecord(x) && typeof x["error"] === "string";
}
function isApiOk(x: unknown): x is ApiOk {
  return (
    isRecord(x) &&
    Array.isArray(x["rows"]) &&
    typeof x["total"] === "number" &&
    typeof x["page"] === "number" &&
    typeof x["pageSize"] === "number"
  );
}

type ColumnKey = keyof RunRow;
type ColumnDef<K extends ColumnKey = ColumnKey> = {
  key: K;
  label: string;
  thClass?: string;
  tdClass?: string;
  numeric?: boolean;
  render?: (row: RunRow) => React.ReactNode;
};

// ---- Column config (labels + per-column classes) ----
const COLUMNS: ColumnDef[] = [
  { key: "job_id",      label: "job_id",      thClass: "text-yellow-300", tdClass: "font-mono text-yellow-700" },
  { key: "finished_at", label: "finished_at" },
  { key: "created_at",  label: "created_at" },
  { key: "dur_minutes", label: "dur_minutes", numeric: true },
  { key: "thr",         label: "thr",         numeric: true },
  { key: "reps",        label: "reps",        numeric: true },
  { key: "max_iter",    label: "max_iter",    numeric: true },
  { key: "d_pi_1",      label: "d_pi_1",      numeric: true },
  { key: "d_pi_2",      label: "d_pi_2",      numeric: true },
  { key: "d_pi_3",      label: "d_pi_3",      numeric: true },
  { key: "d_pi_4",      label: "d_pi_4",      numeric: true },
  { key: "d_m_1",       label: "d_m_1",       numeric: true },
  { key: "d_m_2",       label: "d_m_2",       numeric: true },
  { key: "d_m_3",       label: "d_m_3",       numeric: true },
  { key: "d_m_4",       label: "d_m_4",       numeric: true },
];

const SELECTED_KEY = "emtr:selectedJobId";

export default function RunsPage() {
  const [data, setData] = useState<RunRow[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const [page, setPage] = useState(1);
  const [pageSize, setPageSize] = useState(20);
  const [sort, setSort] = useState("finished_at:desc");
  const [q, setQ] = useState("");
  const [from, setFrom] = useState("");
  const [to, setTo] = useState("");
  const [total, setTotal] = useState(0);

  // New: selected job id
  const [selectedJobId, setSelectedJobId] = useState<string | null>(null);

  // Initialize selection from URL (?job=...) or localStorage
  useEffect(() => {
    if (typeof window === "undefined") return;
    const sp = new URLSearchParams(window.location.search);
    const fromUrl = sp.get("job");
    const fromLS = window.localStorage.getItem(SELECTED_KEY);
    setSelectedJobId(fromUrl ?? fromLS);
  }, []);

  // Persist selection to localStorage + URL + broadcast event
  useEffect(() => {
    if (typeof window === "undefined") return;
    const url = new URL(window.location.href);
    if (selectedJobId && selectedJobId.length > 0) {
      window.localStorage.setItem(SELECTED_KEY, selectedJobId);
      url.searchParams.set("job", selectedJobId);
      // broadcast for other pages/components that care
      const ev: CustomEvent<{ jobId: string }> = new CustomEvent("emtr:selectedJobChanged", {
        detail: { jobId: selectedJobId },
      });
      window.dispatchEvent(ev);
    } else {
      window.localStorage.removeItem(SELECTED_KEY);
      url.searchParams.delete("job");
    }
    window.history.replaceState({}, "", url.toString());
  }, [selectedJobId]);

  async function load() {
    setLoading(true);
    setError(null);
    try {
      const u = new URL("/api/runs", window.location.origin);
      u.searchParams.set("page", String(page));
      u.searchParams.set("pageSize", String(pageSize));
      u.searchParams.set("sort", sort);
      if (q) u.searchParams.set("q", q);
      if (from) u.searchParams.set("from", from);
      if (to) u.searchParams.set("to", to);

      const res = await fetch(u.toString(), { cache: "no-store" });
      const json: unknown = await res.json();

      if (!res.ok) {
        const msg = isApiErr(json) ? json.error : `HTTP ${res.status}`;
        throw new Error(msg);
      }
      if (!isApiOk(json)) throw new Error("Invalid API response");

      setData(json.rows);
      setTotal(json.total);
    } catch (e: unknown) {
      setError(e instanceof Error ? e.message : "Failed to load");
    } finally {
      setLoading(false);
    }
  }

  useEffect(() => {
    void load();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [page, pageSize, sort]);

  function onSearchSubmit(e: React.FormEvent) {
    e.preventDefault();
    setPage(1);
    void load();
  }

  const pages = useMemo(() => Math.max(1, Math.ceil(total / pageSize)), [total, pageSize]);
  const columnCount = 1 + COLUMNS.length; // +1 for the selection column

  const handleSelect = (jobId: string) => setSelectedJobId(jobId);

  return (
    <div className="p-6 max-w-[1400px] mx-auto">
      <h1 className="text-2xl font-bold mb-4">Completed Runs</h1>

      <form onSubmit={onSearchSubmit} className="flex flex-wrap gap-3 items-end mb-4">
        <div className="flex flex-col">
          <label className="text-sm">Search job_id</label>
          <input
            className="border rounded-xl px-3 py-2"
            placeholder="e.g. wasm-17556..."
            value={q}
            onChange={(e) => setQ(e.target.value)}
          />
        </div>
        <div className="flex flex-col">
          <label className="text-sm">From (ISO)</label>
          <input
            type="date"
            className="border rounded-xl px-3 py-2"
            value={from}
            onChange={(e) => setFrom(e.target.value)}
          />
        </div>
        <div className="flex flex-col">
          <label className="text-sm">To (ISO)</label>
          <input
            type="date"
            className="border rounded-xl px-3 py-2"
            value={to}
            onChange={(e) => setTo(e.target.value)}
          />
        </div>
        <div className="flex flex-col">
          <label className="text-sm">Sort</label>
          <select
            className="border rounded-xl px-3 py-2"
            value={sort}
            onChange={(e) => setSort(e.target.value)}
          >
            <option value="finished_at:desc">finished_at ↓</option>
            <option value="finished_at:asc">finished_at ↑</option>
            <option value="created_at:desc">created_at ↓</option>
            <option value="created_at:asc">created_at ↑</option>
            <option value="dur_minutes:desc">duration ↓</option>
            <option value="dur_minutes:asc">duration ↑</option>
            <option value="reps:desc">reps ↓</option>
            <option value="reps:asc">reps ↑</option>
            <option value="thr:desc">thr ↓</option>
            <option value="thr:asc">thr ↑</option>
          </select>
        </div>
        <div className="flex gap-2">
          <button type="submit" className="px-4 py-2 rounded-xl bg-black text-white">
            Apply
          </button>
          <button
            type="button"
            className="px-4 py-2 rounded-xl border"
            onClick={() => {
              setQ("");
              setFrom("");
              setTo("");
              setSort("finished_at:desc");
              setPage(1);
              void load();
            }}
          >
            Reset
          </button>
        </div>
        <div className="ml-auto flex items-center gap-2">
          <label className="text-sm">Page size</label>
          <select
            className="border rounded-xl px-3 py-2"
            value={pageSize}
            onChange={(e) => setPageSize(Number(e.target.value))}
          >
            {[10, 20, 50, 100].map((n) => (
              <option key={n} value={n}>
                {n}
              </option>
            ))}
          </select>
        </div>
      </form>

      <div className="overflow-auto border rounded-2xl">
        <table className="min-w-full text-sm">
          {/* Header */}
          <thead className="bg-black text-white sticky top-0">
            <tr>
              <Th className="text-center w-10">Select</Th>
              {COLUMNS.map((c) => (
                <Th key={String(c.key)} className={c.thClass}>
                  {c.label}
                </Th>
              ))}
            </tr>
          </thead>

          <tbody>
            {loading ? (
              <tr>
                <td colSpan={columnCount} className="p-6 text-center">
                  Loading…
                </td>
              </tr>
            ) : error ? (
              <tr>
                <td colSpan={columnCount} className="p-6 text-center text-red-600">
                  {error}
                </td>
              </tr>
            ) : data.length === 0 ? (
              <tr>
                <td colSpan={columnCount} className="p-6 text-center">
                  No completed runs found
                </td>
              </tr>
            ) : (
              data.map((r) => {
                const selected = r.job_id === selectedJobId;
                return (
                  <tr
                    key={r.job_id}
                    onClick={() => handleSelect(r.job_id)}
                    className={
                      "cursor-pointer " +
                      (selected
                        ? "bg-yellow-50 ring-1 ring-yellow-300 text-black"
                        : "hover:bg-gray-50 hover:text-black")
                    }
                  >
                    {/* Selection radio */}
                    <Td className="text-center align-middle">
                      <input
                        type="radio"
                        name="selectedJob"
                        aria-label={`Select ${r.job_id}`}
                        checked={selected}
                        onChange={() => handleSelect(r.job_id)}
                      />
                    </Td>

                    {/* Data cells */}
                    {COLUMNS.map((c) => {
                      const content = c.render ? c.render(r) : formatCell(r[c.key as keyof RunRow], c.key as keyof RunRow);
                      const numericClass = c.numeric ? "text-right tabular-nums" : "";
                      return (
                        <Td key={String(c.key)} className={`${numericClass} ${c.tdClass ?? ""}`}>
                          {content}
                        </Td>
                      );
                    })}
                  </tr>
                );
              })
            )}
          </tbody>
        </table>
      </div>

      <div className="flex items-center justify-between mt-4">
        <div className="text-sm text-gray-600">
          Page {page} / {pages} · {total} runs
        </div>
        <div className="flex gap-2">
          <button
            className="px-3 py-2 rounded-xl border"
            disabled={page <= 1}
            onClick={() => setPage((p) => Math.max(1, p - 1))}
          >
            Prev
          </button>
          <button
            className="px-3 py-2 rounded-xl border"
            disabled={page >= pages}
            onClick={() => setPage((p) => Math.min(pages, p + 1))}
          >
            Next
          </button>
        </div>
      </div>      
    </div>
  );
}

function Th({ children, className }: { children: React.ReactNode; className?: string }) {
  return <th className={`text-left font-semibold px-3 py-2 whitespace-nowrap ${className ?? ""}`}>{children}</th>;
}

function Td({
  children,
  mono,
  num,
  className,
}: {
  children: React.ReactNode;
  mono?: boolean;
  num?: boolean;
  className?: string;
}) {
  return (
    <td className={`px-3 py-2 ${mono ? "font-mono" : ""} ${num ? "text-right tabular-nums" : ""} ${className ?? ""}`}>
      {children}
    </td>
  );
}

function formatCell(value: RunRow[ColumnKey] | undefined, key: ColumnKey): React.ReactNode {
  if (value == null) return "";
  if (key === "finished_at" || key === "created_at") {
    return fmt(String(value));
  }
  return String(value);
}

function fmt(s: string | null) {
  if (!s) return "";
  const d = new Date(s);
  if (Number.isNaN(d.getTime())) return s;
  return d.toLocaleString();
}
