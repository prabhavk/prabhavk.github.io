// app/mle/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { ProbHists } from "@/components/ProbHists";
import type { PlotlyHTMLElement } from "plotly.js";

type AnyObj = Record<string, unknown>;

/* ---------- new: initialization methods ---------- */
type MethodKey = "dirichlet" | "parsimony" | "hss";
const METHOD_LABEL: Record<MethodKey, string> = {
  dirichlet: "Dirichlet",
  parsimony: "Parsimony",
  hss: "HSS",
};

const SELECTED_JOB_KEY = "emtr:selectedJobId";
const selectedRepKey = (job: string) => `emtr:selectedRep:${job}`;
const selectedRootKey = (job: string) => `emtr:selectedRoot:${job}`;
const selectedMethodKey = (job: string) => `emtr:selectedInitMethod:${job}`;

/* ---------------- number coercion helpers ---------------- */
const asNum = (v: unknown): number | null => {
  if (typeof v === "number" && Number.isFinite(v)) return v;
  if (typeof v === "string" && v.trim() !== "" && Number.isFinite(Number(v))) return Number(v);
  return null;
};
const asNumArray = (x: unknown): number[] | null => {
  if (!Array.isArray(x)) return null;
  const out = x.map(asNum);
  return out.every((v) => v != null) ? (out as number[]) : null;
};
const as4Vec = (x: unknown): number[] | null => {
  const arr = asNumArray(x);
  return arr && arr.length === 4 ? arr : null;
};

const safeParseJSON = (x: unknown): unknown => {
  if (typeof x === "string") {
    try { return JSON.parse(x); } catch { return null; }
  }
  return x;
};

/* --------------- identify single-matrix shapes --------------- */
const looksLikeMatrixObject = (x: unknown): boolean => {
  if (!x || typeof x !== "object") return false;
  const o = x as AnyObj;
  const keys = [
    "rows", "data", "m", "nrow", "ncol", "values", "v", "flat",
    "m11", "m44", "0", "1", "2", "3", "r1", "r4",
  ];
  return keys.some((k) => k in o);
};
const looksLike4x4Array = (x: unknown): boolean =>
  Array.isArray(x) && x.length === 4 && x.every((r) => Array.isArray(r) && (r as unknown[]).length === 4);
const looksLikeFlat16 = (x: unknown): boolean =>
  Array.isArray(x) && x.length === 16 && asNumArray(x) != null;

/* --------------- wrap transitions for ProbHists --------------- */
function toTransMap(val: unknown): Record<string, unknown> | null {
  if (val == null) return null;
  if (typeof val === "object" && !Array.isArray(val)) {
    const o = val as AnyObj;
    const entries = Object.entries(o);
    const anyLooksMatrix = entries.some(
      ([, v]) => looksLike4x4Array(v) || looksLikeFlat16(v) || looksLikeMatrixObject(v)
    );
    if (anyLooksMatrix) return o;
    return { "1": o };
  }
  if (looksLike4x4Array(val) || looksLikeFlat16(val)) return { "1": val as unknown };
  return null;
}

/* ---------------- picking helpers ---------------- */
const pick = (obj: AnyObj, keys: string[]): unknown => {
  for (const k of keys) if (k in obj) return (obj as AnyObj)[k];
  return undefined;
};

/** Validate-and-pick a 4-vector from several candidate keys */
const pick4Vec = (obj: AnyObj, keys: string[]): number[] | null => {
  for (const k of keys) {
    if (!(k in obj)) continue;
    const vec = as4Vec(safeParseJSON((obj as AnyObj)[k]));
    if (vec) return vec;
  }
  return null;
};

/** Validate-and-pick a transitions container from several candidate keys */
const pickTransWrapped = (obj: AnyObj, keys: string[]): Record<string, unknown> | null => {
  for (const k of keys) {
    if (!(k in obj)) continue;
    const wrapped = toTransMap(safeParseJSON((obj as AnyObj)[k]));
    if (wrapped && Object.keys(wrapped).length) return wrapped;
  }
  return null;
};

/* ---------------- API response typing ---------------- */
type ApiMLE = {
  job_id: string;
  method: string;
  rep?: string | number | null;
  reps?: Array<string | number>;
  // canonical names:
  root_prob?: number[] | string | null;
  trans_prob?: unknown;
  // label-only (e.g., "h_5")
  root?: string | null;
  // compatibility fields:
  root_prob_final?: number[] | string | null;
  trans_prob_final?: unknown;
  trans?: unknown;
  D_pi?: number[] | null;
  D_M?: number[] | null;
  error?: string;
};

type ViolinSeriesResp = {
  job_id: string;
  roots: string[];
  root?: string;
};

function okJSON<T>(res: Response): Promise<T> {
  const ct = res.headers.get("content-type") || "";
  if (!ct.toLowerCase().includes("application/json")) {
    return res.text().then((body) => {
      throw new Error(`Expected JSON, got ${ct || "unknown"} (HTTP ${res.status}). ${body.slice(0, 160)}`);
    });
  }
  return res.json() as Promise<T>;
}

/* ---------------- filenames + Plotly helpers ---------------- */
function sanitizeForFilename(s: string): string {
  return s
    .replaceAll(/\s*@\s*/g, "__")
    .replaceAll(/[()]/g, "")
    .replaceAll("|", "_given_")
    .replaceAll(/[^\w\-\.]+/g, "_")
    .replaceAll(/_+/g, "_")
    .replace(/^_+|_+$/g, "");
}
function trimNum(n: number): string {
  if (!Number.isFinite(n)) return "NA";
  const s = n % 1 === 0 ? String(n) : n.toFixed(4);
  return s.replace(/\.?0+$/,"");
}
function paramSlug(name: string, arr?: number[] | null): string {
  if (!arr || arr.length !== 4 || !arr.every(v => Number.isFinite(v))) return `${name}-NA`;
  return `${name}-${arr.map(trimNum).join("_")}`;
}
function guessFigureSuffixFromAnnotationText(t?: string | null): string | null {
  if (!t) return null;
  const stripped = t
    .replace(/\s*\(n=\d+\)\s*$/i, "")
    .replace(/\s*for child node\s+\S+\s*$/i, "");
  const cleaned = sanitizeForFilename(stripped);
  return cleaned ? cleaned.toLowerCase() : null;
}

type GdWithLayout = PlotlyHTMLElement & { layout?: Partial<import("plotly.js").Layout> };
function getFirstAnnotationText(gd: PlotlyHTMLElement): string | null {
  const lay = (gd as GdWithLayout).layout;
  const ann = Array.isArray(lay?.annotations) && lay?.annotations?.length ? lay.annotations[0] : undefined;
  const txt = (ann && typeof (ann as { text?: unknown }).text === "string") ? (ann as { text: string }).text : null;
  return txt ?? null;
}
function dataUrlDownload(filename: string, dataUrl: string): void {
  const a = document.createElement("a");
  a.href = dataUrl;
  a.download = filename;
  a.rel = "noopener";
  document.body.appendChild(a);
  a.click();
  a.remove();
}

/* ------------------------- Page ------------------------- */
export default function MLEinDPage() {
  const [job, setJob] = useState<string>("");
  const [rep, setRep] = useState<string | null>(null);
  const [allReps, setAllReps] = useState<string[]>([]);

  const [roots, setRoots] = useState<string[]>([]);
  const [selRoot, setSelRoot] = useState<string>("");

  /* ---------- new: initialization method state ---------- */
  const [method, setMethod] = useState<MethodKey>("dirichlet");

  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  const [root, setRoot] = useState<number[] | null>(null);                // 4-vector
  const [trans, setTrans] = useState<Record<string, unknown> | null>(null); // map of matrices
  const [rootKey, setRootKey] = useState<string | null>(null);            // label (e.g., "h_5")
  const [alphaPi, setAlphaPi] = useState<number[] | null>(null);
  const [alphaM, setAlphaM] = useState<number[] | null>(null);

  const plotsRef = useRef<HTMLDivElement | null>(null);

  /* ---------- read job + rep + root + method from localStorage ---------- */
  useEffect(() => {
    try {
      const savedJob = localStorage.getItem(SELECTED_JOB_KEY) ?? "";
      setJob(savedJob);
      if (savedJob) {
        const savedRep = localStorage.getItem(selectedRepKey(savedJob));
        setRep(savedRep || null);
        const savedRoot = localStorage.getItem(selectedRootKey(savedJob));
        if (savedRoot) setSelRoot(savedRoot);
        const savedMethod = localStorage.getItem(selectedMethodKey(savedJob)) as MethodKey | null;
        if (savedMethod && (["dirichlet","parsimony","hss"] as MethodKey[]).includes(savedMethod)) {
          setMethod(savedMethod);
        }
      }
    } catch { /* ignore */ }
  }, []);

  /* ---------- persist rep/root/method per job ---------- */
  useEffect(() => {
    if (!job) return;
    try {
      if (rep) localStorage.setItem(selectedRepKey(job), rep);
      else localStorage.removeItem(selectedRepKey(job));
    } catch { /* ignore */ }
  }, [job, rep]);

  useEffect(() => {
    if (!job) return;
    try {
      if (selRoot) localStorage.setItem(selectedRootKey(job), selRoot);
      else localStorage.removeItem(selectedRootKey(job));
    } catch { /* ignore */ }
  }, [job, selRoot]);

  useEffect(() => {
    if (!job) return;
    try {
      if (method) localStorage.setItem(selectedMethodKey(job), method);
      else localStorage.removeItem(selectedMethodKey(job));
    } catch { /* ignore */ }
  }, [job, method]);

  /* ---------- fetch full reps list for this job ---------- */
  const fetchReps = useCallback(async (jobId: string) => {
    try {
      const u = new URL("/api/mle/reps", window.location.origin);
      u.searchParams.set("job_id", jobId);
      const res = await fetch(u.toString(), { cache: "no-store" });
      const j = await res.json().catch(() => null) as { ok?: boolean; reps?: Array<string | number>; error?: string } | null;

      if (res.ok && j && j.ok && Array.isArray(j.reps)) {
        const list = j.reps.map((r) => String(r));
        setAllReps(list);
        if (!rep && list.length) setRep(list[list.length - 1]);
      }
    } catch { /* fallback handled in /api/mle response */ }
  }, [rep]);

  useEffect(() => {
    if (job) void fetchReps(job);
  }, [job, fetchReps]);

  /* ---------- fetch available roots for this job (via /api/violin) ---------- */
  const fetchRoots = useCallback(async (jobId: string) => {
    try {
      const u = new URL("/api/violin", window.location.origin);
      u.searchParams.set("job", jobId);
      const res = await fetch(u.toString(), { cache: "no-store" });
      const json = await okJSON<ViolinSeriesResp | { error: string }>(res);
      if ("error" in json) throw new Error(json.error);

      const list = Array.isArray((json as ViolinSeriesResp).roots) ? (json as ViolinSeriesResp).roots : [];
      setRoots(list);

      // adopt canonical/first root if none selected or invalid
      const canonical = (json as ViolinSeriesResp).root;
      if (!selRoot || (list.length && !list.includes(selRoot))) {
        setSelRoot(canonical ?? list[0] ?? "");
      }
    } catch {
      setRoots([]);
    }
  }, [selRoot]);

  useEffect(() => {
    if (job) void fetchRoots(job);
  }, [job, fetchRoots]);

  /* ---------- load current rep/root/method’s data ---------- */
  const load = useCallback(async () => {
    if (!job) {
      setErr("No job selected. Choose a job on the Precomputed Results page.");
      setRoot(null); setTrans(null); setRootKey(null); setAlphaPi(null); setAlphaM(null);
      setAllReps([]);
      return;
    }

    setLoading(true);
    setErr(null);
    setRoot(null); setTrans(null); setRootKey(null); setAlphaPi(null); setAlphaM(null);

    try {
      const u = new URL("/api/mle", window.location.origin);
      u.searchParams.set("job_id", job);
      u.searchParams.set("method", method);     // <-- use selected initialization method
      if (rep) u.searchParams.set("rep", rep);
      if (selRoot) u.searchParams.set("root", selRoot);

      const res = await fetch(u.toString(), { cache: "no-store" });
      const raw = await res.text();

      let base: unknown = {};
      try { base = raw ? JSON.parse(raw) : {}; } catch {
        throw new Error(`Expected JSON from ${u.pathname}, got invalid body (HTTP ${res.status}). ${raw.slice(0, 180)}`);
      }
      if (!res.ok) throw new Error(`HTTP ${res.status}`);

      const payload = (base as Partial<ApiMLE>) ?? {};
      if (typeof payload.error === "string" && payload.error) throw new Error(payload.error);

      // Reps fallback
      if (!allReps.length && Array.isArray(payload.reps)) {
        const list = payload.reps.map((x) => String(x));
        if (list.length) setAllReps(list);
      }

      const serverRep = payload.rep != null ? String(payload.rep) : null;
      const defaultRep = serverRep || (allReps.length ? allReps[allReps.length - 1] : null);
      if (!rep && defaultRep) { setRep(defaultRep); setLoading(false); return; }

      // Root probabilities (4-vector), if provided by this method
      const rootVec = pick4Vec(payload as AnyObj, [
        "root_prob",          // canonical
        "root_prob_final",    // legacy
        "root_probs",
        "pi",
        "pi_final",
        "root_pi",
      ]);
      setRoot(rootVec);

      // Transitions container (wrapped), if provided
      const transMap = pickTransWrapped(payload as AnyObj, [
        "trans_prob",         // canonical
        "trans_prob_final",   // legacy
        "trans",              // convenience
        "transition",
        "transitions",
        "M_prob_final",
        "P",
        "T",
      ]);
      setTrans(transMap);

      // Root label (server wins; else keep user choice)
    const rKeyServer =
      typeof (payload as AnyObj)["root"] === "string"
      ? ((payload as AnyObj)["root"] as string)
      : (pick(payload as AnyObj, ["rootKey", "root_key"]) as string | undefined);

  setRootKey(selRoot || rKeyServer || null);

      // Dirichlet hyperparameters (likely only when method === 'dirichlet')
      const aPi = as4Vec(payload.D_pi ?? (pick(payload as AnyObj, ["alpha_pi", "dirichlet_pi"]) as unknown)) ?? null;
      const aM  = as4Vec(payload.D_M  ?? (pick(payload as AnyObj, ["alpha_M",  "dirichlet_M"])  as unknown))  ?? null;
      setAlphaPi(aPi);
      setAlphaM(aM);
    } catch (e) {
      setErr(e instanceof Error ? e.message : "Failed to load");
    } finally {
      setLoading(false);
    }
  }, [job, rep, selRoot, method, allReps.length]);

  useEffect(() => { if (job) void load(); }, [job, rep, selRoot, method, load]);

  const rootOK = useMemo(() => Array.isArray(root) && root.length === 4, [root]);
  const transOK = useMemo(() => !!trans && Object.keys(trans).length > 0, [trans]);

  /* ---------------- Download helpers ---------------- */
  function getAllPlots(): PlotlyHTMLElement[] {
    const container = plotsRef.current;
    if (!container) return [];
    return Array.from(container.querySelectorAll<PlotlyHTMLElement>(".js-plotly-plot"));
  }
  function plotFilterRoot(txt: string): boolean {
    return /^p\([acgt]\)/i.test(txt) && !txt.includes("|");
  }
  function plotFilterSingle(txt: string): boolean {
    return /for child node/i.test(txt);
  }
  function plotFilterAgg(txt: string): boolean {
    return /p\(.+\|.+\)/i.test(txt) && /\(n=\d+\)/i.test(txt);
    }

  const dpiSlug = useMemo(() => paramSlug("dpi", alphaPi), [alphaPi]);
  const dmSlug  = useMemo(() => paramSlug("dm",  alphaM),  [alphaM]);

  function prefix(): string {
    const repPart = rep ? `_rep${rep}` : "";
    return `mle_${job || "job"}${repPart}_${method}`;
  }
  function rootNameSlug(): string {
    return sanitizeForFilename(rootKey || "root");
  }

  async function downloadRootPNGs(): Promise<void> {
    const plots = getAllPlots().filter((gd) => {
      const txt = getFirstAnnotationText(gd);
      return !!txt && plotFilterRoot(txt);
    });
    if (!plots.length) { alert("No root-probability plots found."); return; }
    if (!window.Plotly?.toImage) { alert("Plotly export not available."); return; }

    const pre = prefix();
    for (const gd of plots) {
      const txt = getFirstAnnotationText(gd);
      const suffix = guessFigureSuffixFromAnnotationText(txt) || "root_prob";
      const fileName = `${pre}__${suffix}__${dpiSlug}__root-${rootNameSlug()}.png`;
      const url = await window.Plotly.toImage(gd, { format: "png", scale: 2 });
      dataUrlDownload(fileName, url);
    }
  }

  async function downloadSingleEntryPNG(): Promise<void> {
    const gd = getAllPlots().find((g) => {
      const txt = getFirstAnnotationText(g);
      return !!txt && plotFilterSingle(txt);
    });
    if (!gd) { alert("No selected transition-entry plot found."); return; }
    if (!window.Plotly?.toImage) { alert("Plotly export not available."); return; }

    const pre = prefix();
    const txt = getFirstAnnotationText(gd) || "";
    const suffix = guessFigureSuffixFromAnnotationText(txt) || "transition_entry";
    const fileName = `${pre}__${suffix}__${dmSlug}.png`;
    const url = await window.Plotly.toImage(gd, { format: "png", scale: 2 });
    dataUrlDownload(fileName, url);
  }

  async function downloadAggregatedPNGs(): Promise<void> {
    const plots = getAllPlots().filter((gd) => {
      const txt = getFirstAnnotationText(gd);
      return !!txt && plotFilterAgg(txt);
    });
    if (!plots.length) { alert("No aggregated transition-entry plots found."); return; }
    if (!window.Plotly?.toImage) { alert("Plotly export not available."); return; }

    const pre = prefix();
    let idx = 0;
    for (const gd of plots) {
      idx += 1;
      const txt = getFirstAnnotationText(gd);
      const suffix = guessFigureSuffixFromAnnotationText(txt) || `agg_${idx}`;
      const fileName = `${pre}__${suffix}__${dmSlug}.png`;
      const url = await window.Plotly.toImage(gd, { format: "png", scale: 2 });
      dataUrlDownload(fileName, url);
    }
  }

  /* ------------------- UI ------------------- */
  const headerTitle = useMemo(
    () => `Maximum likelihood estimate — ${METHOD_LABEL[method]} initialization`,
    [method]
  );

  return (
    <div className="p-6 max-w-[1400px] mx-auto space-y-4" ref={plotsRef}>
      {/* Header row */}
      <div className="flex flex-wrap items-end gap-3">
        <h1 className="text-2xl text-black font-bold">
          {headerTitle}
        </h1>

        <div className="text-sm ml-2">
          Job: <span className="font-mono">{job || "(none)"} </span>
        </div>

        {/* Repetition selector */}
        <label className="ml-4 text-sm text-black flex items-center gap-2">
          <span>Repetition</span>
          <select
            className="border rounded px-2 py-1 text-sm"
            value={rep ?? ""}
            onChange={(e) => setRep(e.target.value || null)}
            disabled={!allReps.length}
            title={allReps.length ? "Choose repetition" : "No repetitions found"}
          >
            {allReps.length ? (
              allReps.map((r) => (
                <option key={r} value={r}>{r}</option>
              ))
            ) : (
              <option value="">(none)</option>
            )}
          </select>
        </label>

        {/* Root selector */}
        <label className="text-sm text-black flex items-center gap-2">
          <span>Root node</span>
          <select
            className="border rounded px-2 py-1 text-sm"
            value={selRoot}
            onChange={(e) => setSelRoot(e.target.value)}
            disabled={!roots.length}
            title={roots.length ? "Choose root" : "No roots found for this job"}
          >
            {roots.length ? (
              roots.map((n) => (
                <option key={n} value={n}>
                  {n}
                </option>
              ))
            ) : (
              <option value="" disabled>
                {job ? "No roots for this job" : "Select a job first"}
              </option>
            )}
          </select>
        </label>

        {/* NEW: Method selector */}
        <div className="flex items-center gap-2 ml-2">
          <span className="text-sm text-black">Initialization</span>
          <div className="inline-flex rounded-md border overflow-hidden">
            {(["dirichlet","parsimony","hss"] as MethodKey[]).map((m) => (
              <button
                key={m}
                type="button"
                onClick={() => setMethod(m)}
                className={`px-3 py-1.5 text-sm ${
                  method === m ? "bg-black text-white" : "bg-white text-black"
                } ${m !== "hss" ? "border-r" : ""}`}
                title={METHOD_LABEL[m]}
              >
                {METHOD_LABEL[m]}
              </button>
            ))}
          </div>
        </div>

        {rep ? <div className="text-xs text-gray-600">Using rep: <code>{rep}</code></div> : null}

        <button
          type="button"
          onClick={() => void load()}
          className="ml-auto px-4 py-2 rounded bg-gray-600 text-white disabled:opacity-50"
          disabled={!job || loading}
          title="Reload"
        >
          {loading ? "Loading…" : "Reload"}
        </button>
      </div>

      {/* Downloads toolbar */}
      <div className="flex flex-wrap items-center gap-2 rounded p-3 bg-gray-100">
        <div className="text-sm font-medium text-black mr-2">Downloads:</div>
        <button
          type="button"
          onClick={() => void downloadRootPNGs()}
          className="px-3 py-2 rounded border border-black text-black hover:bg-black hover:text-white disabled:opacity-40"
          disabled={!rootOK}
          title="Download the four root-probability plots (PNG)"
        >
          Roots (PNG)
        </button>
        <button
          type="button"
          onClick={() => void downloadSingleEntryPNG()}
          className="px-3 py-2 rounded border border-black text-black hover:bg-black hover:text-white disabled:opacity-40"
          disabled={!transOK}
          title="Download the user-selected transition entry plot (PNG)"
        >
          Selected Entry (PNG)
        </button>
        <button
          type="button"
          onClick={() => void downloadAggregatedPNGs()}
          className="px-3 py-2 rounded border border-black text-black hover:bg-black hover:text-white disabled:opacity-40"
          disabled={!transOK}
          title="Download the 16 aggregated transition-entry plots (PNG)"
        >
          Aggregated (PNG)
        </button>
      </div>

      {loading ? (
        <div className="p-4 border rounded">Loading…</div>
      ) : err ? (
        <div className="p-4 border rounded text-red-600">{err}</div>
      ) : !(rootOK || transOK) ? (
        <div className="p-4 border rounded">
          No {METHOD_LABEL[method]} data available for this job{rep ? ` (rep ${rep})` : ""}.
        </div>
      ) : (
        <div className="space-y-4">
          <ProbHists
            root={rootOK ? root : null}
            trans={transOK ? (trans as Record<string, unknown>) : null}
            rootKey={rootKey ?? undefined}
            rootName={rootKey ?? undefined}
            alphaPi={alphaPi ?? undefined}
            alphaM={alphaM ?? undefined}
            /* Only show Beta CIs when the method is Dirichlet and we have α params */
            showBetaCI={method === "dirichlet" && !!alphaPi}
            betaCILevel={0.95}
            rootHeaderRight={
              <button
                type="button"
                onClick={() => void downloadRootPNGs()}
                className="px-3 py-2 rounded border border-black text-black hover:bg-black hover:text-white disabled:opacity-40"
                disabled={!rootOK}
                title="Download the four root-probability plots (PNG)"
              >
                Download Roots (PNG)
              </button>
            }
            singleHeaderRight={
              <button
                type="button"
                onClick={() => void downloadSingleEntryPNG()}
                className="px-3 py-2 rounded border border-black text-black hover:bg-black hover:text-white disabled:opacity-40"
                disabled={!transOK}
                title="Download the user-selected transition entry plot (PNG)"
              >
                Download Selected Entry (PNG)
              </button>
            }
            aggHeaderRight={
              <button
                type="button"
                onClick={() => void downloadAggregatedPNGs()}
                className="px-3 py-2 rounded border border-black text-black hover:bg-black hover:text-white disabled:opacity-40"
                disabled={!transOK}
                title="Download the 16 aggregated transition-entry plots (PNG)"
              >
                Download Aggregated (PNG)
              </button>
            }
          />
        </div>
      )}
    </div>
  );
}

// window.Plotly type
declare global {
  interface Window {
    Plotly?: {
      toImage: (
        gd: PlotlyHTMLElement,
        opts?: Partial<{ format: "png" | "svg" | "jpeg" | "webp"; width: number; height: number; scale: number }>
      ) => Promise<string>;
    };
  }
}
