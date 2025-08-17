"use client";

import { useSearchParams } from "next/navigation";

export default function ResultsPage() {
  const searchParams = useSearchParams();

  return (
    <main>
      <h1 className="text-2xl font-bold mb-4">EMTR Run Parameters</h1>
      <ul className="list-disc pl-6 text-black">
        <li>Sequence file: {searchParams.get("sequenceFile")}</li>
        <li>Topology file: {searchParams.get("topologyFile")}</li>
        <li>Convergence threshold: {searchParams.get("threshold")}</li>
        <li>Number of repetitions: {searchParams.get("numReps")}</li>
        <li>Max iterations: {searchParams.get("maxIter")}</li>
      </ul>
    </main>
  );
}
