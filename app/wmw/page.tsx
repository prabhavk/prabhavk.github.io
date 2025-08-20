// app/wmw/page.tsx
'use client';

import React, { useEffect, useMemo, useState } from 'react';
import { getSelectedJobId, onSelectedJobChanged } from '@/lib/selectedJob';

type MethodName = 'Parsimony' | 'Dirichlet' | 'SSH';
type PairKey = 'parsimony_vs_dirichlet' | 'parsimony_vs_ssh' | 'dirichlet_vs_ssh';

type PairResult = {
  winner: MethodName | 'none';
  p_value: number; // one-sided p-value for the winning direction
  nA: number;
  nB: number;
};

type ApiResp = {
  job_id: string;
  node?: string;
  alpha: number;
  pairs: Record<PairKey, PairResult>;
};
type ApiErr = { error: string };

function isRecord(x: unknown): x is Record<string, unknown> {
  return typeof x === 'object' && x !== null;
}
function isApiErr(x: unknown): x is ApiErr {
  return isRecord(x) && typeof x['error'] === 'string';
}
function isApiOk(x: unknown): x is ApiResp {
  return isRecord(x) && typeof x['job_id'] === 'string' && isRecord(x['pairs']);
}

// Nodes for the "root" column
const NODES = Array.from({ length: 17 }, (_, i) => `h_${21 + i}`) as readonly string[];


type RowState = {
  node: string;
  pairs: Partial<Record<PairKey, PairResult>>; // filled after "Compute"
};

export default function WmwPage() {
  const [jobId, setJobId] = useState<string | null>(null);
  const [alpha, setAlpha] = useState<number>(0.05);
  const [rows, setRows] = useState<RowState[]>(
    () => NODES.map((node) => ({ node, pairs: {} }))
  );
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Keep rows keyed by node for quick updates
  const rowIndexByNode = useMemo(() => {
    const map = new Map<string, number>();
    rows.forEach((r, i) => map.set(r.node, i));
    return map;
  }, [rows]);

  useEffect(() => {
    setJobId(getSelectedJobId());
    return onSelectedJobChanged((jid) => setJobId(jid));
  }, []);

  function labelWinner(pr?: PairResult): string {
    if (!pr) return '—';
    return pr.winner === 'none' ? '—' : pr.winner;
  }

  function fmtP(pr?: PairResult): string {
    if (!pr) return '—';
    const p = pr.p_value;
    if (!Number.isFinite(p)) return '—';
    if (p < 1e-4) return p.toExponential(2); // e.g., 1.23e-5
    return p.toFixed(4).replace(/0+$/,'').replace(/\.$/,''); // tidy
  }

  async function computeAll() {
    if (!jobId) {
      setError('Select a job in Precomputed Results first.');
      return;
    }
    setError(null);
    setLoading(true);

    try {
      // fetch all nodes in parallel
      const reqs = NODES.map(async (node) => {
        const u = new URL('/api/mwu', window.location.origin);
        u.searchParams.set('job', jobId);
        u.searchParams.set('alpha', String(alpha));
        u.searchParams.set('node', node); // the API will use this to filter emtr_rows by node

        const res = await fetch(u.toString(), { cache: 'no-store' });
        const json: unknown = await res.json();

        if (!res.ok) throw new Error(isApiErr(json) ? json.error : `HTTP ${res.status}`);
        if (!isApiOk(json)) throw new Error('Invalid API response');

        return { node, pairs: json.pairs };
      });

      const results = await Promise.all(reqs);

      // update table state
      setRows((prev) => {
        const next = [...prev];
        for (const r of results) {
          const idx = rowIndexByNode.get(r.node);
          if (idx != null) {
            next[idx] = { node: r.node, pairs: r.pairs };
          }
        }
        return next;
      });
    } catch (e: unknown) {
      setError(e instanceof Error ? e.message : 'Failed to compute');
    } finally {
      setLoading(false);
    }
  }

  return (
    <div className="p-6 max-w-[1100px] mx-auto">
      <h1 className="text-2xl font-bold mb-3">Wilcoxon–Mann–Whitney (per node)</h1>

      <div className="flex flex-wrap items-end gap-3 mb-4">
        <div className="text-sm">
          <div className="text-gray-500">Selected Job</div>
          <div className="font-mono">{jobId ?? '— (select a job in Precomputed Results)'}</div>
        </div>

        <div className="flex items-center gap-2 ml-auto">
          <label className="text-sm">alpha</label>
          <input
            type="number"
            step="0.001"
            min="0.001"
            max="0.5"
            className="border rounded-lg px-3 py-2 w-24"
            value={alpha}
            onChange={(e) => setAlpha(Number(e.target.value))}
          />
          <button
            className="px-4 py-2 rounded-lg bg-black text-white disabled:opacity-50"
            disabled={!jobId || loading}
            onClick={() => void computeAll()}
          >
            {loading ? 'Computing…' : 'Compute'}
          </button>
        </div>
      </div>

      <div className="overflow-auto border rounded-2xl">
        <table className="min-w-full text-sm">
          <thead className="bg-black text-white">
            <tr>
              <th className="text-left font-semibold px-3 py-2">root</th>
              <th className="text-left font-semibold px-3 py-2">Parsimony vs Dirichlet</th>
              <th className="text-left font-semibold px-3 py-2">p</th>
              <th className="text-left font-semibold px-3 py-2">Parsimony vs SSH</th>
              <th className="text-left font-semibold px-3 py-2">p</th>
              <th className="text-left font-semibold px-3 py-2">Dirichlet vs SSH</th>
              <th className="text-left font-semibold px-3 py-2">p</th>
            </tr>
          </thead>
          <tbody>
            {rows.map((r) => (
              <tr key={r.node} className="hover:bg-gray-50">
                <td className="px-3 py-2 font-mono">{r.node}</td>

                <td className="px-3 py-2">{labelWinner(r.pairs.parsimony_vs_dirichlet)}</td>
                <td className="px-3 py-2 tabular-nums">{fmtP(r.pairs.parsimony_vs_dirichlet)}</td>

                <td className="px-3 py-2">{labelWinner(r.pairs.parsimony_vs_ssh)}</td>
                <td className="px-3 py-2 tabular-nums">{fmtP(r.pairs.parsimony_vs_ssh)}</td>

                <td className="px-3 py-2">{labelWinner(r.pairs.dirichlet_vs_ssh)}</td>
                <td className="px-3 py-2 tabular-nums">{fmtP(r.pairs.dirichlet_vs_ssh)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {error && <p className="text-sm text-red-600 mt-3">{error}</p>}

      <p className="text-xs text-gray-500 mt-3">
        Click <em>Compute</em> to populate winners and p-values per node using Wilcoxon–Mann–Whitney tests.
      </p>
    </div>
  );
}
