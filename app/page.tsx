"use client";

import { useState } from "react";
import { useRouter } from "next/navigation";

export default function Home() {
  const router = useRouter();

  const [sequenceFile, setSequenceFile] = useState("");
  const [topologyFile, setTopologyFile] = useState("");
  const [threshold, setThreshold] = useState("");
  const [numReps, setNumReps] = useState("");
  const [maxIter, setMaxIter] = useState("");

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    const params = new URLSearchParams({
      sequenceFile,
      topologyFile,
      threshold,
      numReps,
      maxIter,
    });
    router.push(`/results?${params.toString()}`);
  };

  return (
    <main>
      <h1 className="text-2xl font-bold mb-4">Run EMTR</h1>
      <form onSubmit={handleSubmit} className="space-y-4">
        <input
          type="text"
          placeholder="Sequence file"
          className="border p-2 w-full text-black"
          value={sequenceFile}
          onChange={(e) => setSequenceFile(e.target.value)}
        />
        <input
          type="text"
          placeholder="Topology file"
          className="border p-2 w-full text-black"
          value={topologyFile}
          onChange={(e) => setTopologyFile(e.target.value)}
        />
        <input
          type="text"
          placeholder="Convergence threshold"
          className="border p-2 w-full text-black"
          value={threshold}
          onChange={(e) => setThreshold(e.target.value)}
        />
        <input
          type="number"
          placeholder="Number of repetitions"
          className="border p-2 w-full text-black"
          value={numReps}
          onChange={(e) => setNumReps(e.target.value)}
        />
        <input
          type="number"
          placeholder="Max iterations"
          className="border p-2 w-full text-black"
          value={maxIter}
          onChange={(e) => setMaxIter(e.target.value)}
        />
        <button
          type="submit"
          className="bg-blue-500 text-white px-4 py-2 rounded"
        >
          Submit
        </button>
      </form>
    </main>
  );
}
