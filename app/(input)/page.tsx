// app/(input)/page.tsx — ECDLL (DNA only), 3-stage pattern counts with clamping
"use client";

import React, { useEffect, useRef, useState } from "react";
import { useProgress } from "@/components/ProgressProvider";
import { useProgressStore } from "@/stores/progress";

// ---- Settings schema version (bump when JSON shape changes) ----
const SETTINGS_VERSION = 2; // includes patternCounts

// ---- Validation & guardrails ----
const REQUIRE_ALPHA_GE1_DEFAULT = true;
const ALPHA_MIN = 1e-6;
const ALPHA_MAX = 1e6;
const KAPPA_MAX = 1e8;

const COARSE_MIN = 26;
const COARSE_MAX = 100;
const MEDIUM_MAX = 100;
const FINE_FIXED = 100;

function clamp(n: number, lo: number, hi: number) {
  if (!Number.isFinite(n)) return lo;
  return Math.max(lo, Math.min(hi, Math.trunc(n)));
}

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

type PatternCounts = { coarse: number; medium: number; fine: number };
type LegacyThrPerRegime = { coarse?: number; medium?: number; fine?: number };
type LegacyPatWRegime = { coarse?: number; medium?: number; fine?: number };

type EmtrDNASettings = {
  thr: number;
  reps: number;
  maxIter: number;
  D_pi: number[];
  D_M: number[];
  requireAlphaGe1: boolean;
  includeParsimony?: boolean;
  includeHSS?: boolean;
  patternCounts: PatternCounts;
};

type LegacyEmtrDNASettings = Partial<EmtrDNASettings> & {
  thrPerRegime?: LegacyThrPerRegime;
  patWRegime?: LegacyPatWRegime;
};

type FilesMap = { dna: string };

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

// ---- Worker message types ----
type WorkerInMsg =
  | { __cmd: "setUpsertUiLogging"; enabled: boolean }
  | { __cmd: "setDbLogging"; enabled: boolean }
  | { __cmd: "setTreeUploading"; enabled: boolean }
  | { __cmd: "setApiSecret"; secret: string }
  | {
      __cmd: "runWithJson";
      json: string;
      files: FilesMap;
      dnaSeqBytes?: ArrayBuffer;
      jobId: string;
    };

type WorkerOutMsg =
  | { type: "log"; line: string }
  | { type: "artifact"; name?: string; bytes: number[] }
  | { type: "done"; rc: number };

type TrifleLayer = {
  layer: 0 | 1 | 2;
  iter?: number | null;
  ll_final?: number | null;
  root_prob_final?: [number, number, number, number] | null;
  trans_prob_final?: Record<string, unknown> | null;
  ecd_ll_per_iter?: Record<string, number> | Array<[number, number]> | null;
};

type TrifleBody = {
  job_id: string;
  method: "parsimony" | "dirichlet" | "hss";
  rep: number;
  root: string;
  ll_init?: number | null;
  root_prob_init?: [number, number, number, number] | null;
  trans_prob_init?: Record<string, unknown> | null;
  raw_json?: unknown;
  layers: TrifleLayer[]; // expect 3 entries
};

// Extend the worker message union to include the EM-Trifle message
type ExtendedWorkerOutMsg =
  | WorkerOutMsg
  | { type: "emtrifle"; payload: TrifleBody };

/** Messages that may come from the worker which the page ignores (store handles). */
type WorkerProgressMsg = {
  type: "progress";
  kind:
    | "setup:reps"
    | "setup:edges"
    | "rep"
    | "method:start"
    | "node"
    | "root:done"
    | "layer:done";
  [k: string]: unknown;
};
type WorkerEdgesMsg = { type: "edges"; edges: [string, string, number][] };
type WorkerNewickMsg = { type: "tree_newick"; payload: string };

/** The full set of messages the page might see. */
type KnownWorkerMsg =
  | ExtendedWorkerOutMsg
  | WorkerProgressMsg
  | WorkerEdgesMsg
  | WorkerNewickMsg;

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
  const files: FilesMap = {
    dna: settings.useDefaultSeqs
      ? "/data/RAxML_DNA_test.fa"
      : `/work/${settings.fileNames.dna ?? "dna.fa"}`,
  };

  const payload = {
    job_id: jobId,
    files,
    settings: {
      version: settings.version,
      dna: settings.dna, // includes clamped patternCounts
      flags: settings.flags,
    },
  };

  return { jsonStr: JSON.stringify(payload), files };
}

export default function InputPage() {
  const { start, append } = useProgress();

  const [secret, setSecret] = useState<string>("");

  // Files
  const [DNAsequenceFile, setDNASequenceFile] = useState<File | null>(null);
  const [useDefaultSeqsToggle, setUseDefaultSeqsToggle] = useState<boolean>(false);

  // =========================
  // ECDLL (EM-DNA) controls — single regime
  // =========================
  const [numReps, setNumReps] = useState<string>("30");
  const [maxIter, setMaxIter] = useState<string>("100");
  const [thrDNA, setThrDNA] = useState<string>("0.01");
  const [pi, setPi] = useState<string[]>(["100", "100", "100", "100"]);
  const [M, setM] = useState<string[]>(["100", "2", "2", "2"]);
  const [requireAlphaGe1, setRequireAlphaGe1] = useState<boolean>(REQUIRE_ALPHA_GE1_DEFAULT);

  // Optional pipeline passes
  const [includeParsimony, setIncludeParsimony] = useState<boolean>(true);
  const [includeHSS, setincludeHSS] = useState<boolean>(true);

  // Three-stage site-pattern counts
  const [patCoarse, setPatCoarse] = useState<string>("28");
  const [patMedium, setPatMedium] = useState<string>("49");
  const patFine = String(FINE_FIXED); // fixed

  // Worker flags
  const [enableUiUpsert, setEnableUiUpsert] = useState<boolean>(true);
  const [enableDbLogging, setEnableDbLogging] = useState<boolean>(true);
  const [enableTreeUpload, setEnableTreeUpload] = useState<boolean>(true);

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

  // ---------- helper to push current (validated/clamped) UI inputs into the store ----------
  function pushInputsToStore({
    reps,
    iters,
    thrD,
    D_pi,
    D_M,
    requireAlphaGe1,
    includeParsimony,
    includeHSS,
    pc,
    pm,
    useDefaultSeqs,
    dnaName,
  }: {
    reps: number;
    iters: number;
    thrD: number;
    D_pi: number[];
    D_M: number[];
    requireAlphaGe1: boolean;
    includeParsimony: boolean;
    includeHSS: boolean;
    pc: number;
    pm: number;
    useDefaultSeqs: boolean;
    dnaName: string | null;
  }) {
    // Keep the existing fractions-based wiring (for wedge radii + totalReps)
    useProgressStore.getState().ingestRunInputs({
      reps,
      layerPercGlobal: { middle: pm, top: pc },
    });

    // NEW: store all input params for anyone (e.g., clock) to consume
    useProgressStore.getState().setEmInput({
      thr: thrD,
      reps,
      maxIter: iters,
      D_pi: [D_pi[0], D_pi[1], D_pi[2], D_pi[3]],
      D_M: [D_M[0], D_M[1], D_M[2], D_M[3]],
      requireAlphaGe1,
      includeParsimony,
      includeHSS,
      patternCounts: { coarse: pc, medium: pm, fine: FINE_FIXED },
      useDefaultSeqs,
      dnaFileName: dnaName,
    });
  }

  // ---- Save settings from current UI (without launching job) ----
  function handleSaveSettings() {
    const reps = Number(numReps);
    const iters = Number(maxIter);
    const thrD = Number(thrDNA);
    const D_pi = pi.map(Number);
    const D_M = M.map(Number);

    const errDNA =
      (!Number.isInteger(reps) || reps <= 0 ? "Repetitions must be a positive integer." : null) ||
      (!Number.isInteger(iters) || iters <= 0 ? "Max iterations must be a positive integer." : null) ||
      (!Number.isFinite(thrD) || thrD <= 0 ? "Convergence threshold must be a positive number." : null) ||
      validateDirichlet(D_pi, requireAlphaGe1) ||
      validateDirichlet(D_M, requireAlphaGe1);
    if (errDNA) {
      append(errDNA);
      return;
    }

    // Clamp & enforce ordering: 26 ≤ coarse ≤ 100, coarse ≤ medium ≤ 100
    const pc = clamp(Number(patCoarse), COARSE_MIN, COARSE_MAX);
    const pm = clamp(Number(patMedium), pc, MEDIUM_MAX);
    if (pc !== Number(patCoarse)) setPatCoarse(String(pc));
    if (pm !== Number(patMedium)) setPatMedium(String(pm));

    const defaultDNA = "/data/RAxML_DNA_test.fa";
    const useDefaultSeqs = useDefaultSeqsToggle || !DNAsequenceFile;
    const dnaName = !useDefaultSeqs ? DNAsequenceFile?.name ?? null : defaultDNA;

    // NEW: push to store
    pushInputsToStore({
      reps, iters, thrD, D_pi, D_M,
      requireAlphaGe1, includeParsimony, includeHSS,
      pc, pm, useDefaultSeqs, dnaName
    });

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
        patternCounts: { coarse: pc, medium: pm, fine: FINE_FIXED },
      },
      useDefaultSeqs,
      fileNames: { dna: dnaName },
      flags: { enableUiUpsert, enableDbLogging, enableTreeUpload },
    };

    downloadSettingsJSON(settings);
    append("Saved settings JSON.");
  }

  // ---- Load settings JSON into UI ----
  async function handleLoadSettings(file: File) {
    try {
      const text = await file.text();
      const parsed: unknown = JSON.parse(text);
      const s: EmtrSettingsPartial = parsed as EmtrSettingsPartial;

      if (!s || typeof s !== "object") throw new Error("Invalid JSON.");
      if (isNumber(s.version) && s.version > SETTINGS_VERSION) {
        append(`Warning: settings version ${s.version} is newer than supported (${SETTINGS_VERSION}).`);
      }

      let reps: number | undefined;
      let iters: number | undefined;
      let thrD: number | undefined;
      let D_pi: number[] | undefined;
      let D_M: number[] | undefined;
      let reqGe1: boolean | undefined;
      let incParsimony: boolean | undefined;
      let incHSS: boolean | undefined;
      let pc: number | undefined;
      let pm: number | undefined;

      if (s.dna) {
        const d = s.dna;

        if (isNumber(d.maxIter)) { setMaxIter(String(d.maxIter)); iters = d.maxIter; }
        if (isNumber(d.reps))    { setNumReps(String(d.reps));    reps = d.reps; }

        let thrCandidate: number | undefined = undefined;
        if (isNumber(d.thr)) thrCandidate = d.thr;
        else if (d.thrPerRegime) {
          const tr = d.thrPerRegime;
          if (isNumber(tr.fine))   thrCandidate = tr.fine;
          else if (isNumber(tr.medium)) thrCandidate = tr.medium;
          else if (isNumber(tr.coarse)) thrCandidate = tr.coarse;
        }
        if (isNumber(thrCandidate)) { setThrDNA(String(thrCandidate)); thrD = thrCandidate; }

        if (Array.isArray(d.D_pi) && d.D_pi.length === 4 && d.D_pi.every(isNumber)) {
          setPi(d.D_pi.map((x) => String(x)));
          D_pi = d.D_pi;
        }
        if (Array.isArray(d.D_M) && d.D_M.length === 4 && d.D_M.every(isNumber)) {
          setM(d.D_M.map((x) => String(x)));
          D_M = d.D_M;
        }

        if (typeof d.requireAlphaGe1 === "boolean") { setRequireAlphaGe1(d.requireAlphaGe1); reqGe1 = d.requireAlphaGe1; }
        if (typeof d.includeParsimony === "boolean") { setIncludeParsimony(d.includeParsimony); incParsimony = d.includeParsimony; }
        if (typeof d.includeHSS === "boolean") { setincludeHSS(d.includeHSS); incHSS = d.includeHSS; }

        // Clamp loaded pattern counts with new rules
        const pcIn = d.patternCounts?.coarse ?? d.patWRegime?.coarse ?? COARSE_MIN;
        const pcClamped = clamp(Number(pcIn), COARSE_MIN, COARSE_MAX);
        setPatCoarse(String(pcClamped));
        pc = pcClamped;

        const pmIn = d.patternCounts?.medium ?? d.patWRegime?.medium ?? pcClamped;
        const pmClamped = clamp(Number(pmIn), pcClamped, MEDIUM_MAX);
        setPatMedium(String(pmClamped));
        pm = pmClamped;
        // fine fixed
      }

      if (typeof s.useDefaultSeqs === "boolean") setUseDefaultSeqsToggle(s.useDefaultSeqs);

      if (s.flags) {
        if (typeof s.flags.enableUiUpsert === "boolean") setEnableUiUpsert(s.flags.enableUiUpsert);
        if (typeof s.flags.enableDbLogging === "boolean") setEnableDbLogging(s.flags.enableDbLogging);
        if (typeof s.flags.enableTreeUpload === "boolean") setEnableTreeUpload(s.flags.enableTreeUpload);
      }

      // NEW: push to store (use current UI after applying loaded values)
      const repsFinal = Number(reps ?? numReps);
      const itersFinal = Number(iters ?? maxIter);
      const thrDFinal = Number(thrD ?? thrDNA);
      const D_piFinal = (D_pi ?? pi.map(Number));
      const D_MFinal  = (D_M ?? M.map(Number));
      const reqGe1Final = (reqGe1 ?? requireAlphaGe1);
      const incParFinal = (incParsimony ?? includeParsimony);
      const incHSSFinal = (incHSS ?? includeHSS);
      const pcFinal = clamp(Number(pc ?? patCoarse), COARSE_MIN, COARSE_MAX);
      const pmFinal = clamp(Number(pm ?? patMedium), pcFinal, MEDIUM_MAX);

      const defaultDNA = "/data/RAxML_DNA_test.fa";
      const useDefaultSeqs = useDefaultSeqsToggle || !DNAsequenceFile;
      const dnaName = !useDefaultSeqs ? DNAsequenceFile?.name ?? null : defaultDNA;

      pushInputsToStore({
        reps: repsFinal,
        iters: itersFinal,
        thrD: thrDFinal,
        D_pi: D_piFinal,
        D_M: D_MFinal,
        requireAlphaGe1: reqGe1Final,
        includeParsimony: incParFinal,
        includeHSS: incHSSFinal,
        pc: pcFinal,
        pm: pmFinal,
        useDefaultSeqs,
        dnaName,
      });

      append(`Loaded settings from ${file.name}. Review and click "start".`);
    } catch (err: unknown) {
      const message = err instanceof Error ? err.message : String(err);
      append(`Failed to load settings: ${message}`);
    }
  }

  // ---- Submit / run ----
  async function handleSubmit(e: React.FormEvent<HTMLFormElement>) {
    e.preventDefault();

    const reps = Number(numReps);
    const iters = Number(maxIter);
    const thrD = Number(thrDNA);

    if (!Number.isInteger(reps) || reps <= 0) {
      append("Repetitions must be a positive integer.");
      return;
    }
    if (!Number.isInteger(iters) || iters <= 0) {
      append("Max iterations must be a positive integer.");
      return;
    }
    if (!Number.isFinite(thrD) || thrD <= 0) {
      append("Convergence threshold must be a positive number.");
      return;
    }

    // Enforce new constraints on submit as well
    const pcRaw = Number(patCoarse);
    const pc = clamp(pcRaw, COARSE_MIN, COARSE_MAX);

    const pmRaw = Number(patMedium);
    const pm = clamp(pmRaw, pc, MEDIUM_MAX);

    if (pc !== pcRaw) setPatCoarse(String(pc));
    if (pm !== pmRaw) setPatMedium(String(pm));

    // fractions used by the clock for wedge radii + totalReps
    useProgressStore.getState().ingestRunInputs({
      reps,
      layerPercGlobal: { middle: pm, top: pc },
    });

    const D_pi = pi.map(Number);
    const D_M = M.map(Number);
    {
      const err = validateDirichlet(D_pi, requireAlphaGe1) || validateDirichlet(D_M, requireAlphaGe1);
      if (err) {
        append(`${err}`);
        return;
      }
    }

    // NEW: also persist all inputs into the store (so the clock & others can read them)
    {
      const defaultDNA = "/data/RAxML_DNA_test.fa";
      const useDefaultSeqs = useDefaultSeqsToggle || !DNAsequenceFile;
      const dnaName = !useDefaultSeqs ? DNAsequenceFile?.name ?? null : defaultDNA;

      useProgressStore.getState().setEmInput({
        thr: thrD,
        reps,
        maxIter: iters,
        D_pi: [D_pi[0], D_pi[1], D_pi[2], D_pi[3]],
        D_M: [D_M[0], D_M[1], D_M[2], D_M[3]],
        requireAlphaGe1,
        includeParsimony,
        includeHSS,
        patternCounts: { coarse: pc, medium: pm, fine: FINE_FIXED },
        useDefaultSeqs,
        dnaFileName: dnaName,
      });
    }

    // Create the job and get job_id
    let jobId: string;
    try {
      const createRes = await fetch("/api/jobs", {
        method: "POST",
        headers: {
          "content-type": "application/json",
          ...(secret ? { "x-emtr-secret": secret } : {}),
        },
        body: JSON.stringify({
          thr: thrD,
          reps,
          max_iter: iters,
          D_pi,
          D_M,
          pattern_counts: { coarse: pc, medium: pm, fine: FINE_FIXED }, // persisted for convenience
        }),
      });
      if (!createRes.ok) {
        const text = await createRes.text();
        throw new Error(`Create job failed: ${createRes.status} ${text}`);
      }
      const data: { job_id: string } = await createRes.json();
      jobId = data.job_id;
    } catch (err: unknown) {
      append(`Failed to create job: ${err instanceof Error ? err.message : String(err)}`);
      return;
    }

    try {
      localStorage.setItem("emtr:selectedJobId", jobId);
    } catch {}

    start(jobId);
    append(`jobId = ${jobId}`);

    const defaultDNA = "/data/RAxML_DNA_test.fa";
    const useDefaultSeqs = useDefaultSeqsToggle || !DNAsequenceFile;

    const parts: string[] = [];
    if (!useDefaultSeqs && DNAsequenceFile) parts.push(`DNA=${DNAsequenceFile.name}`);
    append(
      `starting run... ${parts.join(", ") || `no files → default: DNA=${defaultDNA}`}; ` +
        `EMreps=${reps}, maxEMiter=${iters}; thr=${thrD}; ` +
        `patterns (coarse,medium,fine)=(${pc},${pm},${FINE_FIXED})`
    );

    const v = Date.now();
    const w = new Worker(`/wasm/worker.js?v=${v}`, { type: "module" });
    useProgressStore.getState().attachWorker(w);
    workerRef.current = w;

    w.onmessage = async (ev: MessageEvent<KnownWorkerMsg>) => {
      const msg = ev.data;

      switch (msg.type) {
        case "log": {
          append(msg.line);
          return;
        }
        case "artifact": {
          const blob = new Blob([new Uint8Array(msg.bytes)], { type: "application/octet-stream" });
          const url = URL.createObjectURL(blob);
          const a = document.createElement("a");
          const name = msg.name ?? "emtr_output.bin";
          a.href = url; a.download = name; a.click();
          URL.revokeObjectURL(url);
          return;
        }
        case "emtrifle": {
          try {
            const t = msg.payload;
            const layers =
              Array.isArray(t.layers)
                ? t.layers.map(L => `L${L.layer}:iter=${L.iter ?? "-"} ll=${L.ll_final ?? "-"}`).join(", ")
                : "(no layers)";
            append(`[EM_Trifle] job=${t.job_id} method=${t.method} root=${t.root} rep=${t.rep} | ${layers}`);
            if (localStorage.getItem("emtr.debug.trifle") === "1") {
              append(`[EM_Trifle JSON] ${JSON.stringify(t)}`);
            }
          } catch {}
          return;
        }
        case "done": {
          append(msg.rc === 0 ? "Finished" : "Error");
          w.terminate();
          workerRef.current = null;
          return;
        }

        // These message types are consumed by the global store via attachWorker(w).
        // We intentionally ignore them here to avoid double-processing.
        case "progress":
        case "edges":
        case "tree_newick":
        default:
          return;
      }
    };

    const post = (message: WorkerInMsg) => w.postMessage(message);

    post({ __cmd: "setUpsertUiLogging", enabled: enableUiUpsert });
    post({ __cmd: "setDbLogging", enabled: enableDbLogging });
    post({ __cmd: "setTreeUploading", enabled: enableTreeUpload });
    post({ __cmd: "setApiSecret", secret });

    const dnaBuf: ArrayBuffer | undefined =
      !useDefaultSeqs && DNAsequenceFile ? await DNAsequenceFile.arrayBuffer() : undefined;

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
        patternCounts: { coarse: pc, medium: pm, fine: FINE_FIXED },
      },
      useDefaultSeqs,
      fileNames: {
        dna: !useDefaultSeqs ? DNAsequenceFile?.name ?? null : defaultDNA,
      },
      flags: { enableUiUpsert, enableDbLogging, enableTreeUpload },
    };

    const { jsonStr, files } = buildRunJSON(jobId, settings);

    post({
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
            onChange={(e: React.ChangeEvent<HTMLInputElement>) =>
              setDNASequenceFile(e.target.files?.[0] ?? null)
            }
            className="block w-full border p-2 text-black"
            disabled={useDefaultSeqsToggle}
          />
        </div>

        {/* Settings */}
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
                onChange={(e: React.ChangeEvent<HTMLInputElement>) => setNumReps(e.target.value)}
                placeholder="e.g., 30"
                min={1}
                step={1}
              />
            </div>
            <div>
              <label className="block text-sm font-medium text-black">Most iterations per layer</label>
              <input
                type="number"
                className="w-full border p-2 text-black"
                value={maxIter}
                onChange={(e: React.ChangeEvent<HTMLInputElement>) => setMaxIter(e.target.value)}
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
              onChange={(e: React.ChangeEvent<HTMLInputElement>) => setThrDNA(e.target.value)}
              placeholder="e.g., 0.01"
              min={0}
            />
          </div>

          {/* Three-stage layer counts with ordering & clamping */}
          <div className="space-y-2 pt-2">
            <label className="block text-sm font-medium text-black">
              Cumulative fractional weight of site patterns per layer
            </label>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-2">
              {/* Coarse */}
              <div>
                <div className="text-xs text-black/70 pb-1">Coarse</div>
                <input
                  type="number"
                  min={COARSE_MIN}
                  max={COARSE_MAX}
                  step={1}
                  className="w-full border p-2 text-black"
                  value={patCoarse}
                  onChange={(e: React.ChangeEvent<HTMLInputElement>) => {
                    const newCoarse = clamp(Number(e.target.value), COARSE_MIN, COARSE_MAX);
                    setPatCoarse(String(newCoarse));
                    // ensure medium >= coarse
                    const mediumNum = Number(patMedium);
                    if (!Number.isFinite(mediumNum) || mediumNum < newCoarse) {
                      setPatMedium(String(newCoarse));
                    }
                  }}
                  placeholder="26"
                />
                <p className="text-[11px] text-white mt-1">Allowed: {COARSE_MIN}–{COARSE_MAX}</p>
              </div>

              {/* Medium */}
              <div>
                <div className="text-xs text-black/70 pb-1">Medium</div>
                <input
                  type="number"
                  min={Number(patCoarse) || COARSE_MIN}
                  max={MEDIUM_MAX}
                  step={1}
                  className="w-full border p-2 text-black"
                  value={patMedium}
                  onChange={(e: React.ChangeEvent<HTMLInputElement>) => {
                    const coarseNum = clamp(Number(patCoarse), COARSE_MIN, COARSE_MAX);
                    const v = clamp(Number(e.target.value), coarseNum, MEDIUM_MAX);
                    setPatMedium(String(v));
                  }}
                  placeholder="89"
                />
                <p className="text-[11px] text-white mt-1">Allowed: {patCoarse || COARSE_MIN}–{MEDIUM_MAX}</p>
              </div>

              {/* Fine */}
              <div>
                <div className="text-xs text-black/70 pb-1">Fine</div>
                <input
                  type="number"
                  className="w-full border p-2 text-black bg-gray-100"
                  value={patFine}
                  disabled
                />
                <p className="text-[11px] text-white mt-1">Fixed: {FINE_FIXED}</p>
              </div>
            </div>
          </div>

          {/* Dirichlet α for π */}
          <div className="space-y-2">
            <label className="block text-sm font-medium text-black">Dirichlet α for π</label>
            <div className="grid grid-cols-4 gap-2">
              {pi.map((v, i) => (
                <input
                  key={i}
                  type="number"
                  min={0}
                  step="any"
                  className="w-full border p-2 text-black"
                  value={v}
                  onChange={(e: React.ChangeEvent<HTMLInputElement>) => {
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
                  min={0}
                  step="any"
                  className="w-full border p-2 text-black"
                  value={v}
                  onChange={(e: React.ChangeEvent<HTMLInputElement>) => handleMChange(i, e.target.value)}
                  title={i === 0 ? "α1 (independent)" : "α2–α4 are linked; editing one updates the others"}
                />
              ))}
            </div>
          </div>
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
              onChange={(e: React.ChangeEvent<HTMLInputElement>) => {
                const file = e.target.files?.[0];
                if (file) handleLoadSettings(file);
                e.currentTarget.value = "";
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
              onChange={(e: React.ChangeEvent<HTMLInputElement>) => setSecret(e.target.value)}
              placeholder=""
              autoComplete="current-password"
            />
          </label>

          <button type="submit" className="px-4 py-2 rounded-md bg-yellow-700 text-yellow-300 hover:bg-yellow-800">
            start
          </button>
        </div>
      </form>
    </main>
  );
}
