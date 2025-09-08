// app/(input)/page.tsx — ECDLL (DNA only), single-regime; save/load/run wired; API Secret passthrough
"use client";

import React, { useEffect, useRef, useState } from "react";
import { useProgress } from "@/components/ProgressProvider";

// ---- Settings schema version (bump when JSON shape changes) ----
const SETTINGS_VERSION = 1;

// ---- Validation & guardrails ----
const REQUIRE_ALPHA_GE1_DEFAULT = true;
const ALPHA_MIN = 1e-6;
const ALPHA_MAX = 1e6;
const KAPPA_MAX = 1e8;

function validateDirichlet(a: number[], requireAlphaGe1: boolean): string | null {
  if (!Array.isArray(a) || a.length !== 4) return "Must have exactly 4 α values.";
  for (const x of a) {
    if (!Number.isFinite(x)) return "α values must be finite numbers.";
    if (x <= 0) return "All α must be > 0.";
    if (x < ALPHA_MIN) return `α too small (< ${ALPHA_MIN}) — may be numerically unstable.`;
    if (x > ALPHA_MAX) return `α too large (> ${ALPHA_MAX}) — prior-like weight may dominate.`;
    if (requireAlphaGe1 && x < 1) return "Use α ≥ 1 to avoid boundary/degenerate estimates.";
  }
  const kappa = a[0] + a[1] + a[2] + a[3];
  if (kappa > KAPPA_MAX) return `Sum of α (κ) is too large (> ${KAPPA_MAX}).`;
  return null;
}

type LegacyThrPerRegime = { coarse?: number; medium?: number; fine?: number };
type LegacyPatWRegime  = { coarse?: number; medium?: number; fine?: number };

type EmtrDNASettings = {
  // ECDLL (DNA), single regime (one global threshold)
  thr: number;
  reps: number;
  maxIter: number;
  D_pi: number[]; // length 4
  D_M: number[];  // length 4, with UI tying of 2–4
  requireAlphaGe1: boolean;
  includeParsimony?: boolean;
  includeHSS?: boolean;
};

type LegacyEmtrDNASettings = Partial<EmtrDNASettings> & {
  thrPerRegime?: LegacyThrPerRegime; // legacy only
  patWRegime?: LegacyPatWRegime;     // legacy only
};

type EmtrSettings = {
  version: number;
  dna: EmtrDNASettings;
  useDefaultSeqs: boolean;
  fileNames: { dna: string | null };
  flags: {
    enableUiUpsert: boolean;
    enableDbLogging: boolean;
    enableTreeUpload: boolean;
  };
};

type EmtrSettingsPartial = {
  version?: number;
  dna?: LegacyEmtrDNASettings;
  useDefaultSeqs?: boolean;
  fileNames?: { dna?: string | null };
  flags?: Partial<EmtrSettings["flags"]>;
};

function isNumber(x: unknown): x is number {
  return typeof x === "number" && Number.isFinite(x);
}

// ---- Helpers ----
function downloadSettingsJSON(settings: EmtrSettings, filename = "emtr_settings.json") {
  const blob = new Blob([JSON.stringify(settings, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function buildRunJSON(jobId: string, settings: EmtrSettings) {
  // Choose *target* MEMFS paths the worker/WASM will write to
  const files = {
    dna: settings.useDefaultSeqs
      ? "/data/RAxML_DNA_test.fa"
      : `/work/${settings.fileNames.dna ?? "dna.fa"}`,
  };

  const payload = {
    job_id: jobId,
    files,
    settings: {
      version: settings.version,
      dna: settings.dna,
      flags: settings.flags,
    },
  };

  return { jsonStr: JSON.stringify(payload), files };
}

export default function InputPage() {
  const { start, append } = useProgress();

  // Single secret used for API auth (forwarded to worker → header)
  const [secret, setSecret] = useState("");

  // Files
  const [DNAsequenceFile, setDNASequenceFile] = useState<File | null>(null);
  const [useDefaultSeqsToggle, setUseDefaultSeqsToggle] = useState(false);

  // =========================
  // ECDLL (EM-DNA) controls — single regime
  // =========================
  const [numReps, setNumReps] = useState("10");
  const [maxIter, setMaxIter] = useState("100");
  const [thrDNA, setThrDNA]   = useState("0.01");
  const [pi, setPi]           = useState(["100", "100", "100", "100"]);
  const [M, setM]             = useState(["100", "2", "2", "2"]);
  const [requireAlphaGe1, setRequireAlphaGe1] = useState(REQUIRE_ALPHA_GE1_DEFAULT);

  // Optional pipeline passes
  const [includeParsimony, setIncludeParsimony] = useState(true);
  const [includeHSS, setincludeHSS] = useState(true);

  // Worker flags
  const [enableUiUpsert, setEnableUiUpsert] = useState(true);
  const [enableDbLogging, setEnableDbLogging] = useState(true);
  const [enableTreeUpload, setEnableTreeUpload] = useState(true);

  function handleMChange(index: number, value: string) {
    setM((prev) => {
      const next = [...prev];
      if (index === 0) next[0] = value;
      else next[1] = next[2] = next[3] = value; // lock α2–α4
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

  // ---- Save settings from current UI (without launching job) ----
  function handleSaveSettings() {
    const reps = Number(numReps);
    const iters = Number(maxIter);
    const thrD = Number(thrDNA);
    const D_pi = pi.map(Number);
    const D_M = M.map(Number);

    // validations
    const errDNA =
      (!Number.isInteger(reps) || reps <= 0 ? "Repetitions must be a positive integer." : null) ||
      (!Number.isInteger(iters) || iters <= 0 ? "Max iterations must be a positive integer." : null) ||
      (!Number.isFinite(thrD) || thrD <= 0 ? "Convergence threshold must be a positive number." : null) ||
      validateDirichlet(D_pi, requireAlphaGe1) ||
      validateDirichlet(D_M, requireAlphaGe1);
    if (errDNA) { append(errDNA); return; }

    const defaultDNA = "/data/RAxML_DNA_test.fa";
    const useDefaultSeqs = useDefaultSeqsToggle || !DNAsequenceFile;

    const settings: EmtrSettings = {
      version: SETTINGS_VERSION,
      dna: {
        thr: thrD,
        reps,
        maxIter: iters,
        D_pi,
        D_M,
        requireAlphaGe1,
        includeParsimony,
        includeHSS,
      },
      useDefaultSeqs,
      fileNames: {
        dna: !useDefaultSeqs ? DNAsequenceFile?.name ?? null : defaultDNA,
      },
      flags: { enableUiUpsert, enableDbLogging, enableTreeUpload },
    };

    downloadSettingsJSON(settings);
    append("Saved settings JSON.");
  }

  // ---- Load settings JSON into UI ----
  async function handleLoadSettings(file: File) {
    try {
      const text = await file.text();
      const parsed = JSON.parse(text) as unknown;

      // Narrow to our partial shape (no `any`)
      const s = parsed as EmtrSettingsPartial;

      if (!s || typeof s !== "object") throw new Error("Invalid JSON.");

      if (isNumber(s.version) && s.version > SETTINGS_VERSION) {
        append(`Warning: settings version ${s.version} is newer than supported (${SETTINGS_VERSION}).`);
      }

      if (s.dna) {
        const d = s.dna;

        if (isNumber(d.maxIter)) setMaxIter(String(d.maxIter));
        if (isNumber(d.reps))    setNumReps(String(d.reps));

        // single/global thr (fallback to legacy per-regime values if present)
        let thrCandidate: number | undefined = undefined;
        if (isNumber(d.thr)) thrCandidate = d.thr;
        else if (d.thrPerRegime) {
          const tr = d.thrPerRegime;
          if (isNumber(tr.fine))   thrCandidate = tr.fine;
          else if (isNumber(tr.medium)) thrCandidate = tr.medium;
          else if (isNumber(tr.coarse)) thrCandidate = tr.coarse;
        }
        if (isNumber(thrCandidate)) setThrDNA(String(thrCandidate));

        if (Array.isArray(d.D_pi) && d.D_pi.length === 4 && d.D_pi.every(isNumber)) {
          setPi(d.D_pi.map((x) => String(x)));
        }
        if (Array.isArray(d.D_M) && d.D_M.length === 4 && d.D_M.every(isNumber)) {
          setM(d.D_M.map((x) => String(x)));
        }

        if (typeof d.requireAlphaGe1 === "boolean") setRequireAlphaGe1(d.requireAlphaGe1);
        if (typeof d.includeParsimony === "boolean") setIncludeParsimony(d.includeParsimony);
        if (typeof d.includeHSS === "boolean") setincludeHSS(d.includeHSS);
      }

      if (typeof s.useDefaultSeqs === "boolean") setUseDefaultSeqsToggle(s.useDefaultSeqs);

      if (s.flags) {
        if (typeof s.flags.enableUiUpsert === "boolean") setEnableUiUpsert(s.flags.enableUiUpsert);
        if (typeof s.flags.enableDbLogging === "boolean") setEnableDbLogging(s.flags.enableDbLogging);
        if (typeof s.flags.enableTreeUpload === "boolean") setEnableTreeUpload(s.flags.enableTreeUpload);
      }

      append(`Loaded settings from ${file.name}. Review and click "start qrep".`);
    } catch (err: unknown) {
      const message = err instanceof Error ? err.message : String(err);
      append(`Failed to load settings: ${message}`);
    }
  }

  // ---- Submit / run ----
  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();

    const reps  = Number(numReps);
    const iters = Number(maxIter);
    const thrD  = Number(thrDNA);

    // DNA/ECDLL validation
    if (!Number.isInteger(reps) || reps <= 0) { append("Repetitions must be a positive integer."); return; }
    if (!Number.isInteger(iters) || iters <= 0) { append("Max iterations must be a positive integer."); return; }
    if (!Number.isFinite(thrD) || thrD <= 0) { append("Convergence threshold must be a positive number."); return; }

    const D_pi = pi.map(Number);
    const D_M  = M.map(Number);
    {
      const err = validateDirichlet(D_pi, requireAlphaGe1) || validateDirichlet(D_M, requireAlphaGe1);
      if (err) { append(`${err}`); return; }
    }
// Create the job on the server and get the single source-of-truth job_id
    let jobId: string;
    try {
      const createRes = await fetch("/api/jobs", {
        method: "POST",
        headers: {
          "content-type": "application/json",
          ...(secret ? { "x-emtr-secret": secret } : {}),
        },
        body: JSON.stringify({
          // store core params so emtr_jobs has them immediately
          thr: thrD,
          reps,
          max_iter: iters,
          D_pi,
          D_M,
        }),
      });
      if (!createRes.ok) {
        const text = await createRes.text();
        throw new Error(`Create job failed: ${createRes.status} ${text}`);
      }
      const data = await createRes.json();
      jobId = data.job_id;
    } catch (err: unknown) {
      append(`Failed to create job: ${err instanceof Error ? err.message : String(err)}`);
      return;
    }

    try { localStorage.setItem("emtr:selectedJobId", jobId); } catch {}
    start(jobId);
    append(`jobId = ${jobId}`);

    const defaultDNA = "/data/RAxML_DNA_test.fa";
    const useDefaultSeqs = useDefaultSeqsToggle || !DNAsequenceFile;

    const parts: string[] = [];
    if (!useDefaultSeqs && DNAsequenceFile) parts.push(`DNA=${DNAsequenceFile.name}`);
    append(
      `starting run... ${parts.join(", ") || `no files → default: DNA=${defaultDNA}`}; ` +
      `EMreps=${reps}, maxEMiter=${iters}; ` +
      `thr=${thrD}`
    );

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
        const name = msg.name || "emtr_output.bin";
        a.href = url; a.download = name; a.click();
        URL.revokeObjectURL(url);
      }
      if (msg.type === "done") {
        append(msg.rc === 0 ? "Finished" : "Error");
        w.terminate();
        workerRef.current = null;
      }
    };

    // Worker flags
    w.postMessage({ __cmd: "setUpsertUiLogging", enabled: enableUiUpsert });
    w.postMessage({ __cmd: "setDbLogging",       enabled: enableDbLogging });
    w.postMessage({ __cmd: "setTreeUploading",   enabled: enableTreeUpload });

    // Forward secret to worker so it can add the header
    w.postMessage({ __cmd: "setApiSecret", secret });

    // Prepare sequence buffer (only if user uploaded)
    const dnaBuf = !useDefaultSeqs && DNAsequenceFile ? await DNAsequenceFile.arrayBuffer() : undefined;

    // Construct grouped settings object to send to worker
    const settings: EmtrSettings = {
      version: SETTINGS_VERSION,
      dna: {
        thr: thrD,
        reps,
        maxIter: iters,
        D_pi,
        D_M,
        requireAlphaGe1,
        includeParsimony,
        includeHSS,
      },
      useDefaultSeqs,
      fileNames: {
        dna: !useDefaultSeqs ? DNAsequenceFile?.name ?? null : defaultDNA,
      },
      flags: { enableUiUpsert, enableDbLogging, enableTreeUpload },
    };

    // (Optional) persist a copy per run
    try { downloadSettingsJSON(settings, `emtr_settings_${jobId}.json`); } catch {}

    // Build the single JSON string for C++
    const { jsonStr, files } = buildRunJSON(jobId, settings);

    // Post *once* to the worker: JSON string + raw bytes (if any)
    w.postMessage({
    __cmd: "runWithJson",
      json: jsonStr,
      files,
      dnaSeqBytes: dnaBuf,
      jobId, 
    });
  }

  return (
    <main>
      <h1 className="text-2xl text-black font-bold mb-4">Input</h1>

      <form onSubmit={handleSubmit} className="space-y-6">
        {/* DNA file */}
        <div className="space-y-2">
          <label className="block text-sm font-medium text-black">DNA sequence file (.fa)</label>
          <input
            type="file"
            accept=".fa,.fasta,.phy,.phylip,.txt"
            onChange={(e) => setDNASequenceFile(e.target.files?.[0] ?? null)}
            className="block w-full border p-2 text-black"
            disabled={useDefaultSeqsToggle}
          />
        </div>

        {/* Use default dataset toggle */}
        <div className="flex items-center gap-3">
          <input
            id="use-defaults"
            type="checkbox"
            className="h-4 w-4"
            checked={useDefaultSeqsToggle}
            onChange={(e) => setUseDefaultSeqsToggle(e.target.checked)}
          />
          <label htmlFor="use-defaults" className="text-sm text-black">
            Use default dataset path (ignore upload)
          </label>
        </div>

        {/* ===================== Settings for EM-DNA (ECDLL) ===================== */}
        <div className="rounded-2xl border border-stone-300 p-4 space-y-4 bg-white">
          <h2 className="text-lg font-semibold text-black">Settings for EMDNA</h2>

          {/* Numeric parameters (reps & maxIter) */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <div>
              <label className="block text-sm font-medium text-black">Repetitions</label>
              <input
                type="number"
                className="w-full border p-2 text-black"
                value={numReps}
                onChange={(e) => setNumReps(e.target.value)}
                placeholder="e.g., 10"
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
                placeholder="e.g., 100"
                min={1}
                step={1}
              />
            </div>
          </div>

          {/* Single convergence threshold */}
          <div>
            <label className="block text-sm font-medium text-black">Convergence threshold</label>
            <input
              type="number"
              step="any"
              className="w-full border p-2 text-black mt-1"
              value={thrDNA}
              onChange={(e)=>setThrDNA(e.target.value)}
              placeholder="e.g., 0.01"
              min={0}
            />
          </div>

          {/* Stability toggle (α ≥ 1) */}
          <div className="flex items-center gap-3">
            <input
              id="alpha-ge1"
              type="checkbox"
              className="h-4 w-4"
              checked={requireAlphaGe1}
              onChange={(e) => setRequireAlphaGe1(e.target.checked)}
            />
            <label htmlFor="alpha-ge1" className="text-sm text-black">
              Enforce α ≥ 1 (stability)
            </label>
          </div>

          {/* Dirichlet α for π */}
          <div className="space-y-2">
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
          </div>

          {/* Dirichlet α for rows of P */}
          <div className="space-y-2">
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
              <strong>Constraint:</strong> α₂, α₃ and α₄ are tied
            </p>
          </div>

          {/* Include Parsimony / HSS */}
          <div className="flex items-center gap-3 pt-1">
            <button
              type="button"
              aria-pressed={includeParsimony}
              onClick={() => setIncludeParsimony(v => !v)}
              className={
                "px-3 py-1.5 rounded border text-sm " +
                (includeParsimony
                  ? "bg-[#FF6B3D] text-white border-gray-800"
                  : "bg-white text-black border-gray-400 hover:bg-gray-50")
              }
              title="Toggle inclusion of a parsimony pass before EM"
            >
              Include parsimony
            </button>

            <button
              type="button"
              onClick={() => setincludeHSS(!includeHSS)}
              className={
                "px-3 py-1.5 rounded border text-sm " +
                (includeHSS
                  ? "bg-[#BB1E10] text-white border-gray-700"
                  : "bg-white text-black border-gray-400 hover:bg-gray-50")
              }
              title="Toggle inclusion of HSS step"
            >
              Include HSS
            </button>
          </div>
        </div>

        {/* Worker flags */}
        <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
          <label className="flex items-center gap-2 text-black">
            <input type="checkbox" className="h-4 w-4" checked={enableUiUpsert} onChange={(e) => setEnableUiUpsert(e.target.checked)} />
            Upsert UI logs
          </label>
          <label className="flex items-center gap-2 text-black">
            <input type="checkbox" className="h-4 w-4" checked={enableDbLogging} onChange={(e) => setEnableDbLogging(e.target.checked)} />
            DB logging
          </label>
          <label className="flex items-center gap-2 text-black">
            <input type="checkbox" className="h-4 w-4" checked={enableTreeUpload} onChange={(e) => setEnableTreeUpload(e.target.checked)} />
            Upload trees
          </label>
        </div>

        {/* Settings I/O */}
        <div className="flex flex-wrap items-center gap-3">
          <button type="button" onClick={handleSaveSettings} className="border px-3 py-1 rounded text-black">
            Save settings (.json)
          </button>

          <label className="cursor-pointer border px-3 py-1 rounded text-black">
            Load settings (.json)
            <input
              type="file"
              accept="application/json,.json"
              className="hidden"
              onChange={(e) => {
                const file = e.target.files?.[0];
                if (file) handleLoadSettings(file);
                (e.target as HTMLInputElement).value = ""; // allow reloading same file later
              }}
            />
          </label>
        </div>

        {/* API secret field + CTA */}
        <div className="flex items-center gap-3">
          <label className="text-sm text-black flex items-center gap-2">
            <span className="whitespace-nowrap">enter secret</span>
            <input
              type="password"
              className="border rounded px-2 py-1 text-black"
              value={secret}
              onChange={(e) => setSecret(e.target.value)}
              placeholder=""
              autoComplete="current-password"
            />
          </label>

          <button
            type="submit"
            className="px-4 py-2 rounded-md border bg-white text-black hover:bg-gray-50 focus:outline-none focus-visible:ring-2 focus-visible:ring-zinc-400"
          >
            start
          </button>
        </div>
      </form>
    </main>
  );
}
