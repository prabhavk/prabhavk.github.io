// app/(input)/page.tsx
"use client";

import React, { useEffect, useRef, useState } from "react";
import { useProgress } from "@/components/ProgressProvider";

const REQUIRE_MAP_SAFE = true;
const ALPHA_MIN = 1e-6;
const ALPHA_MAX = 1e6;
const KAPPA_MAX = 1e8;

function validateDirichlet(a: number[], requireMapSafe = REQUIRE_MAP_SAFE): string | null {
  if (!Array.isArray(a) || a.length !== 4) return "Must have exactly 4 α values.";
  for (const x of a) {
    if (!Number.isFinite(x)) return "α values must be finite numbers.";
    if (x <= 0) return "All α must be > 0.";
    if (x < ALPHA_MIN) return `α too small (< ${ALPHA_MIN}) — may be numerically unstable.`;
    if (x > ALPHA_MAX) return `α too large (> ${ALPHA_MAX}) — prior may dominate.`;
    if (requireMapSafe && x < 1) return "For MAP updates, use α ≥ 1 to avoid degenerate numerators.";
  }
  const kappa = a[0] + a[1] + a[2] + a[3];
  if (kappa > KAPPA_MAX) return `Sum of α (κ) is too large (> ${KAPPA_MAX}).`;
  return null;
}

export default function InputPage() {
  const { start, append } = useProgress();

  const [DNAsequenceFile, setDNASequenceFile] = useState<File | null>(null);
  const [AAsequenceFile,  setAASequenceFile]  = useState<File | null>(null);

  const [threshold, setThreshold] = useState("0.001");
  const [numReps,   setNumReps]   = useState("25");
  const [maxIter,   setMaxIter]   = useState("1000");
  const [pi, setPi] = useState(["100", "100", "100", "100"]);
  const [M,  setM]  = useState(["100", "2", "2", "2"]); // α2–α4 linked

  function handleMChange(index: number, value: string) {
    setM((prev) => {
      const next = [...prev];
      if (index === 0) next[0] = value;
      else next[1] = next[2] = next[3] = value;
      return next;
    });
  }

  const workerRef = useRef<Worker | null>(null);

  useEffect(() => {
    return () => {
      workerRef.current?.terminate();
      workerRef.current = null;
    };
  }, []);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();

    // Require at least one of DNA/AA
    if (!DNAsequenceFile && !AAsequenceFile) {
      append("Please select at least one sequence file (DNA and/or AA).");
      return;
    }

    const thr   = Number(threshold);
    const reps  = Number(numReps);
    const iters = Number(maxIter);

    if (!Number.isFinite(thr) || thr <= 0) {
      append("Convergence threshold must be a positive number.");
      return;
    }
    if (!Number.isInteger(reps) || reps <= 0) {
      append("Repetitions must be a positive integer.");
      return;
    }
    if (!Number.isInteger(iters) || iters <= 0) {
      append("Max iterations must be a positive integer.");
      return;
    }

    const D_pi = pi.map(Number);
    const D_M  = M.map(Number);
    const err = validateDirichlet(D_pi) || validateDirichlet(D_M);
    if (err) { append(`${err}`); return; }

    const jobId = `wasm-${Date.now()}`;
    try { localStorage.setItem("emtr:selectedJobId", jobId); } catch {}

    start(jobId);
    append(`jobId = ${jobId}`);

    // Friendly summary line
    const parts: string[] = [];
    if (DNAsequenceFile) parts.push(`DNA=${DNAsequenceFile.name}`);
    if (AAsequenceFile)  parts.push(`AA=${AAsequenceFile.name}`);
    append(`pouring batter... ${parts.join(", ") || "no files"}; EMthresh=${thr}, EMreps=${reps}, maxEMiter=${iters}`);

    const v = Date.now();
    const w = new Worker(`/wasm/worker.js?v=${v}`, { type: "module" });
    workerRef.current = w;

    w.onmessage = (ev) => {
      const msg = ev.data;
      if (msg.type === "log") append(msg.line);
      if (msg.type === "artifact") {
        const blob = new Blob([new Uint8Array(msg.bytes)], { type: "application/octet-stream" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = msg.name || "emtr_output.bin";
        a.click();
        URL.revokeObjectURL(url);
      }
      if (msg.type === "done") {
        append(msg.rc === 0 ? "Finished" : "Error");
        w.terminate();
        workerRef.current = null;
      }
    };

    // (Dev) verbose worker logging
    w.postMessage({ __cmd: "setUpsertUiLogging", enabled: true });
    w.postMessage({ __cmd: "setDbLogging", enabled: true });
    w.postMessage({ __cmd: "setTreeUploading", enabled: true });

    // Always send both buffers if present
    const dnaBuf = DNAsequenceFile ? await DNAsequenceFile.arrayBuffer() : undefined;
    const aaBuf  = AAsequenceFile  ? await AAsequenceFile.arrayBuffer()  : undefined;

    w.postMessage({
      params: { thr, reps, maxIter: iters, D_pi, D_M },
      dnaSeqBytes: dnaBuf,
      aaSeqBytes:  aaBuf,
      jobId,
      fileNames: {
        dna: DNAsequenceFile?.name || null,
        aa:  AAsequenceFile?.name  || null,
      },
    });
  }

  return (
    <main>
      <h1 className="text-2xl text-black font-bold mb-4">Input</h1>

      <form onSubmit={handleSubmit} className="space-y-4">
        {/* DNA file */}
        <div className="space-y-2">
          <label className="block text-sm font-medium text-black">DNA sequence file (.fa)</label>
          <input
            type="file"
            accept=".fa,.fasta,.phy,.phylip,.txt"
            onChange={(e) => setDNASequenceFile(e.target.files?.[0] ?? null)}
            className="block w-full border p-2 text-black"
          />
        </div>

        {/* AA file */}
        <div className="space-y-2">
          <label className="block text-sm font-medium text-black">AA sequence file (.fa)</label>
          <input
            type="file"
            accept=".faa,.fa,.fasta,.txt"
            onChange={(e) => setAASequenceFile(e.target.files?.[0] ?? null)}
            className="block w-full border p-2 text-black"
          />
        </div>

        {/* Numeric parameters */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
          <div>
            <label className="block text-sm font-medium text-black">Convergence threshold</label>
            <input
              type="number"
              step="any"
              className="w-full border p-2 text-black"
              value={threshold}
              onChange={(e) => setThreshold(e.target.value)}
              placeholder="e.g., 0.01"
            />
          </div>
          <div>
            <label className="block text-sm font-medium text-black">Repetitions</label>
            <input
              type="number"
              className="w-full border p-2 text-black"
              value={numReps}
              onChange={(e) => setNumReps(e.target.value)}
              placeholder="e.g., 50"
              min={1}
              step={1}
            />
          </div>
          <div>
            <label className="block text-sm font-medium text-black">Max iterations</label>
            <input
              type="number"
              className="w-full border p-2 text-black"
              value={maxIter}
              onChange={(e) => setMaxIter(e.target.value)}
              placeholder="e.g., 1000"
              min={1}
              step={1}
            />
          </div>

          {/* Dirichlet priors (π) */}
          <div className="space-y-2 md:col-span-3">
            <label className="block text-sm font-medium text-black">Dirichlet α for π</label>
            <div className="grid grid-cols-4 gap-2">
              {pi.map((v, i) => (
                <input
                  key={i}
                  type="number"
                  min="0"
                  step="any"
                  className="w-full border p-2 text-black"
                  value={v}
                  onChange={(e) => {
                    const next = [...pi];
                    next[i] = e.target.value;
                    setPi(next);
                  }}
                />
              ))}
            </div>
            <p className="text-xs text-gray-300"> </p>
          </div>

          {/* Dirichlet priors (M) with α₂=α₃=α₄ lock */}
          <div className="space-y-2 md:col-span-3">
            <label className="block text-sm font-medium text-black">Dirichlet α for rows of P</label>
            <div className="grid grid-cols-4 gap-2">
              {M.map((v, i) => (
                <input
                  key={i}
                  type="number"
                  min="0"
                  step="any"
                  className="w-full border p-2 text-black"
                  value={v}
                  onChange={(e) => handleMChange(i, e.target.value)}
                  title={i === 0 ? "α1 (independent)" : "α2–α4 are linked; editing one updates the others"}
                />
              ))}
            </div>
            <p className="text-xs text-black">
              <strong>Constraint:</strong> α₂, α₃, α₄ are tied—changing one updates the others.
            </p>
          </div>
        </div>

        <button className="bg-yellow-800 text-black px-4 py-2 rounded hover:bg-yellow-800 hover:text-yellow-300">
          start qrep
        </button>
      </form>
    </main>
  );
}
