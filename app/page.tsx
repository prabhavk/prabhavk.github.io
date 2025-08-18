"use client";

import { useEffect, useRef, useState } from "react";
import { useProgress } from "../components/ProgressProvider";

type Method = "main" | "dirichlet" | "parsimony" | "ssh";

export default function InputPage() {
  const { start, append } = useProgress();

  // real files instead of strings
  const [sequenceFile, setSequenceFile] = useState<File | null>(null);
  const [topologyFile, setTopologyFile] = useState<File | null>(null);

  // params
  const [threshold, setThreshold] = useState("0.001");
  const [numReps, setNumReps] = useState("50");
  const [maxIter, setMaxIter] = useState("200");
  const [seqFormat, setSeqFormat] = useState("phylip");
  const [method, setMethod] = useState<Method>("main");

  const workerRef = useRef<Worker | null>(null);

  // clean up worker when navigating away
  useEffect(() => {
    return () => {
      workerRef.current?.terminate();
      workerRef.current = null;
    };
  }, []);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();

    if (!sequenceFile || !topologyFile) {
      append("❗ Please select both the sequence and topology files.");
      return;
    }
    const thr = Number(threshold);
    const reps = Number(numReps);
    const iters = Number(maxIter);
    if (!Number.isFinite(thr) || thr <= 0) {
      append("❗ Convergence threshold must be a positive number.");
      return;
    }
    if (!Number.isInteger(reps) || reps <= 0) {
      append("❗ Repetitions must be a positive integer.");
      return;
    }
    if (!Number.isInteger(iters) || iters <= 0) {
      append("❗ Max iterations must be a positive integer.");
      return;
    }

    const jobId = `wasm-${Date.now()}`;
    start(jobId);
    append(`🆔 jobId = ${jobId}`);
    append(
      `🌲 Starting EMTR (WASM): seq=${sequenceFile.name}, topo=${topologyFile.name}, thr=${thr}, reps=${reps}, maxIter=${iters}, format=${seqFormat}, method=${method}`
    );

    // spin up worker (cache-busted to avoid stale worker/module in dev)
    const v = Date.now();
    const w = new Worker(`/wasm/worker.js?v=${v}`, { type: "module" });
    workerRef.current = w;

    w.onmessage = (ev) => {
      const msg = ev.data;
      if (msg.type === "log") append(msg.line);
      if (msg.type === "artifact") {
        // receive an output file from WASM & trigger download
        const blob = new Blob([new Uint8Array(msg.bytes)], {
          type: "application/octet-stream",
        });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = msg.name || "emtr_output.bin";
        a.click();
        URL.revokeObjectURL(url);
      }
      if (msg.type === "done") {
        append(msg.rc === 0 ? "🎄 Finished" : "❌ Error");
        w.terminate();
        workerRef.current = null;
      }
    };

    // send inputs + params
    w.postMessage({
      params: {
        method,
        seqFormat,
        thr,
        reps,
        maxIter: iters,
      },
      seqBytes: await sequenceFile.arrayBuffer(),
      topoBytes: await topologyFile.arrayBuffer(),
      jobId,
    });
  }

  return (
    <main>
      <h1 className="text-2xl font-bold mb-4">Input</h1>

      <form onSubmit={handleSubmit} className="space-y-4">
        {/* Files */}
        <div className="space-y-2">
          <label className="block text-sm font-medium text-white">
            Sequence file (.phyx / phylip)
          </label>
          <input
            type="file"
            onChange={(e) => setSequenceFile(e.target.files?.[0] ?? null)}
            className="block w-full border p-2 text-white"
          />
        </div>

        <div className="space-y-2">
          <label className="block text-sm font-medium text-white">
            Topology file (.csv)
          </label>
          <input
            type="file"
            onChange={(e) => setTopologyFile(e.target.files?.[0] ?? null)}
            className="block w-full border p-2 text-white"
          />
        </div>

        {/* Parameters */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
          <div>
            <label className="block text-sm font-medium text-white">
              Convergence threshold
            </label>
            <input
              className="w-full border p-2 text-white"
              value={threshold}
              onChange={(e) => setThreshold(e.target.value)}
              placeholder="e.g., 0.001"
              inputMode="decimal"
            />
          </div>
          <div>
            <label className="block text-sm font-medium text-white">
              Repetitions
            </label>
            <input
              type="number"
              className="w-full border p-2 text-white"
              value={numReps}
              onChange={(e) => setNumReps(e.target.value)}
              placeholder="e.g., 50"
              min={1}
            />
          </div>
          <div>
            <label className="block text-sm font-medium text-white">
              Max iterations
            </label>
            <input
              type="number"
              className="w-full border p-2 text-white"
              value={maxIter}
              onChange={(e) => setMaxIter(e.target.value)}
              placeholder="e.g., 200"
              min={1}
            />
          </div>
        </div>

        {/* Format & Method */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
          <div>
            <label className="block text-sm font-medium text-white">
              Sequence format
            </label>
            <select
              className="w-full border p-2 text-white"
              value={seqFormat}
              onChange={(e) => setSeqFormat(e.target.value)}
            >
              <option value="phylip">phylip</option>
              {/* <option value="fasta" disabled>fasta (not implemented)</option> */}
            </select>
          </div>
          <div>
            <label className="block text-sm font-medium text-white">
              Method
            </label>
            <select
              className="w-full border p-2 text-white"
              value={method}
              onChange={(e) => setMethod(e.target.value as Method)}
            >
              <option value="main">Main</option>
              <option value="dirichlet">Dirichlet</option>
              <option value="parsimony">Parsimony</option>
              <option value="ssh">SSH</option>
            </select>
          </div>
        </div>

        <button className="bg-gray-600 text-white px-4 py-2 rounded">
          Start EMTR 
        </button>
      </form>
    </main>
  );
}
