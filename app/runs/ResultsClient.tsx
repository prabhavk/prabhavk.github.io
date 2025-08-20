"use client";

import { useSearchParams } from "next/navigation";

export default function ResultsClient() {
  const sp = useSearchParams();
  // read values you passed via the query, e.g.
  const seq = sp.get("sequence") ?? "";
  const topo = sp.get("topology") ?? "";
  const thr = sp.get("thr") ?? "";
  const reps = sp.get("reps") ?? "";
  const maxIter = sp.get("maxIter") ?? "";

  return (
    <div className="space-y-2">
      <h1 className="text-2xl font-bold">Results</h1>
      <div className="text-sm">
        <div>Sequence: {seq}</div>
        <div>Topology: {topo}</div>
        <div>Threshold: {thr}</div>
        <div>Repetitions: {reps}</div>
        <div>Max Iter: {maxIter}</div>
      </div>
    </div>
  );
}
