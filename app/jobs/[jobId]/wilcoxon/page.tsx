"use client";

import { useEffect, useState } from "react";
import { useParams } from "next/navigation";

type PairRes = { U: number; p: number; q?: number; n1: number; n2: number };
type RowRes = {
  root: string;
  parsimony_vs_dirichlet: PairRes;
  parsimony_vs_ssh: PairRes;
  dirichlet_vs_ssh: PairRes;
};

export default function WilcoxonPage() {
  const { jobId } = useParams<{ jobId: string }>();
  const [rows, setRows] = useState<RowRes[]>([]);
  const [err, setErr] = useState<string | null>(null);

  useEffect(() => {
    if (!jobId) return;
    fetch(`/api/jobs/${encodeURIComponent(jobId)}/wilcoxon`)
      .then((r) => r.json())
      .then((j) => {
        if (!j.ok) setErr(j.error || "Error");
        else setRows(j.results || []);
      })
      .catch((e) => setErr(String(e)));
  }, [jobId]);

  if (!jobId) return <div className="p-4">Loading…</div>;
  if (err) return <div className="p-4 text-red-400">Error: {err}</div>;

  return (
    <main className="p-4">
      <h1 className="text-xl font-bold mb-4">Wilcoxon rank-sum (job {jobId})</h1>
      <div className="overflow-x-auto">
        <table className="min-w-full text-sm border">
          <thead>
            <tr className="bg-gray-800">
              <th className="p-2 text-left">Root</th>
              <th className="p-2 text-left">Parsimony vs Dirichlet</th>
              <th className="p-2 text-left">Parsimony vs SSH</th>
              <th className="p-2 text-left">Dirichlet vs SSH</th>
            </tr>
          </thead>
          <tbody>
            {rows.map((r) => (
              <tr key={r.root} className="border-t">
                <td className="p-2 font-mono">{r.root}</td>
                {(["parsimony_vs_dirichlet","parsimony_vs_ssh","dirichlet_vs_ssh"] as const).map((k) => {
                  const v = r[k];
                  const fmt = (x?: number) => (Number.isFinite(x!) ? x!.toExponential(3) : "NA");
                  return (
                    <td key={k} className="p-2">
                      <div>U={fmt(v.U)} | p={fmt(v.p)}{v.q!=null ? ` | q=${fmt(v.q)}` : ""}</div>
                      <div className="text-gray-400">n=({v.n1},{v.n2})</div>
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </main>
  );
}
