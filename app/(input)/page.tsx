"use client";

import { useEffect, useRef, useState } from "react";
import { useProgress } from "@/components/ProgressProvider";

const REQUIRE_MAP_SAFE = true; // set to false if your EM uses posterior means, not MAP
const ALPHA_MIN = 1e-6;       // numeric guardrail
const ALPHA_MAX = 1e6;        // numeric guardrail
const KAPPA_MAX = 1e8;        // sum(α) guardrail

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

  const [sequenceFile, setSequenceFile] = useState<File | null>(null);
  const [topologyFile, setTopologyFile] = useState<File | null>(null);

  const [threshold, setThreshold] = useState("0.01");
  const [numReps, setNumReps] = useState("50");
  const [maxIter, setMaxIter] = useState("1000");

  // Dirichlet priors
  const [pi, setPi] = useState(["100", "100", "100", "100"]); // α for π
  const [M, setM] = useState(["100", "2", "2", "2"]);         // α for one row of M

  const workerRef = useRef<Worker | null>(null);

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

    // Dirichlet α validation
    const D_pi = pi.map(Number);
    const D_M  = M.map(Number);
    let err = validateDirichlet(D_pi);
    if (!err) err = validateDirichlet(D_M);
    if (err) {
      append(`❗ ${err}`);
      return;
    }

    const jobId = `wasm-${Date.now()}`;
    start(jobId);
    append(`🆔 jobId = ${jobId}`);
    append(
      `🌲 Starting EMTR (WASM): seq=${sequenceFile.name}, topo=${topologyFile.name}, thr=${thr}, reps=${reps}, maxIter=${iters}`
    );

    const v = Date.now();
    const w = new Worker(`/wasm/worker.js?v=${v}`, { type: "module" });
    workerRef.current = w;

    w.onmessage = (ev) => {
      const msg = ev.data;
      if (msg.type === "log") append(msg.line);
      if (msg.type === "artifact") {
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

    // Send params (incl. Dirichlet α’s)
    w.postMessage({
      params: {
        thr,
        reps,
        maxIter: iters,
        D_pi, // [pi1, pi2, pi3, pi4]
        D_M,  // [M1, M2, M3, M4]
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
            Sequence file (.phyx)
          </label>
          <input
            type="file"
            accept=".phyx"
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
            accept=".csv"
            onChange={(e) => setTopologyFile(e.target.files?.[0] ?? null)}
            className="block w-full border p-2 text-white"
          />
        </div>

        {/* Numeric parameters */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
          <div>
            <label className="block text-sm font-medium text-white">
              Convergence threshold
            </label>
            <input
              className="w-full border p-2 text-white"
              value={threshold}
              onChange={(e) => setThreshold(e.target.value)}
              placeholder="e.g., 0.01"
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

          {/* Dirichlet priors (span full width on md+) */}
          <div className="space-y-2 md:col-span-3">
            <label className="block text-sm font-medium text-white">
              Dirichlet α for π
            </label>
            <div className="grid grid-cols-4 gap-2">
              {pi.map((v, i) => (
                <input
                  key={i}
                  type="number"
                  min="0"
                  step="any"
                  className="w-full border p-2 text-white"
                  value={v}
                  onChange={(e) => {
                    const next = [...pi];
                    next[i] = e.target.value;
                    setPi(next);
                  }}
                />
              ))}
            </div>
            <p className="text-xs text-gray-300">
              Use α ≥ 1 (defaults: 100,100,100,100).
            </p>
          </div>

          <div className="space-y-2 md:col-span-3">
            <label className="block text-sm font-medium text-white">
              Dirichlet α for rows of M
            </label>
            <div className="grid grid-cols-4 gap-2">
              {M.map((v, i) => (
                <input
                  key={i}
                  type="number"
                  min="0"
                  step="any"
                  className="w-full border p-2 text-white"
                  value={v}
                  onChange={(e) => {
                    const next = [...M];
                    next[i] = e.target.value;
                    setM(next);
                  }}
                />
              ))}
            </div>
            <p className="text-xs text-gray-300">
              Defaults: 100,2,2,2 (strong self-transition prior).
            </p>
          </div>
        </div>

        <button className="bg-gray-400 text-white px-4 py-2 rounded  hover:bg-gray-500 hover:text-white">
          Start EMTR
        </button>
      </form>
    </main>
  );
}
