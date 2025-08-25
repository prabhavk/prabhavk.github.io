// app/mle/page.tsx
"use client";

import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { ProbHists } from "@/components/ProbHists";
import type { PlotlyHTMLElement } from "plotly.js";

type AnyObj = Record<string, unknown>;

const SELECTED_JOB_KEY = "emtr:selectedJobId";
const selectedRepKey = (job: string) => `emtr:selectedRep:${job}`;

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
  Array.isArray(x) &&
  x.length === 4 &&
  x.every((r) => Array.isArray(r) && (r as unknown[]).length === 4);
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
    if (anyLooksMatrix) return o;      // already a map of matrices
    return { "1": o };                 // single matrix-like object -> wrap
  }

  if (looksLike4x4Array(val) || looksLikeFlat16(val)) return { "1": val as unknown };
  return null;
}

/* ---------------- picking helpers ---------------- */
const pick = (obj: AnyObj, keys: string[]): unknown => {
  for (const k of keys) if (k in obj) return obj[k];
  return undefined;
};

/* ---------------- API response typing ---------------- */
type ApiMLE = {
  job_id: string;
  method: string;
  rep?: string | number | null;
  reps?: Array<string | number>;
  root?: unknown;
  trans?: unknown;
  root_name?: string | null;
  D_pi?: number[] | null;
  D_M?: number[] | null;
  error?: string;
};
function isApiMLE(x: unknown): x is ApiMLE {
  return !!x && typeof x === "object" && "job_id" in (x as Record<string, unknown>);
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
/** Drop trailing "(n=...)" and "... for child node <id>" bits from annotation text. */
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

  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  const [root, setRoot] = useState<number[] | null>(null);
  const [trans, setTrans] = useState<Record<string, unknown> | null>(null);
  const [rootKey, setRootKey] = useState<string | null>(null);
  const [alphaPi, setAlphaPi] = useState<number[] | null>(null);
  const [alphaM, setAlphaM] = useState<number[] | null>(null);

  // container around all plots for the querySelector lookup
  const plotsRef = useRef<HTMLDivElement | null>(null);

  // read job + rep from localStorage
  useEffect(() => {
    try {
      const savedJob = localStorage.getItem(SELECTED_JOB_KEY) ?? "";
      setJob(savedJob);
      if (savedJob) {
        const savedRep = localStorage.getItem(selectedRepKey(savedJob));
        setRep(savedRep || null);
      }
    } catch { /* ignore */ }
  }, []);

  // persist rep per job
  useEffect(() => {
    if (!job) return;
    try {
      if (rep) localStorage.setItem(selectedRepKey(job), rep);
      else localStorage.removeItem(selectedRepKey(job));
    } catch { /* ignore */ }
  }, [job, rep]);

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
      u.searchParams.set("method", "dirichlet");
      if (rep) u.searchParams.set("rep", rep);

      const res = await fetch(u.toString(), { cache: "no-store" });
      const raw = await res.text();

      let base: unknown = {};
      try { base = raw ? JSON.parse(raw) : {}; } catch {
        throw new Error(`Expected JSON from ${u.pathname}, got invalid body (HTTP ${res.status}). ${raw.slice(0, 180)}`);
      }
      if (!res.ok) throw new Error(`HTTP ${res.status}`);

      const payload = (base as Partial<ApiMLE>) ?? {};
      if (typeof payload.error === "string" && payload.error) throw new Error(payload.error);

      const serverReps = Array.isArray(payload.reps) ? payload.reps.map((x) => String(x)) : [];
      setAllReps(serverReps);

      const serverRep = payload.rep != null ? String(payload.rep) : null;
      const defaultRep = serverRep || (serverReps.length ? serverReps[serverReps.length - 1] : null);
      if (!rep && defaultRep) { setRep(defaultRep); setLoading(false); return; }

      const rootCandidate = payload.root ?? pick(payload as AnyObj, [
        "root_prob_final", "root_prob", "root_probs", "pi", "pi_final", "root_pi",
      ]);
      setRoot(as4Vec(safeParseJSON(rootCandidate)));

      const transCandidate = payload.trans ?? pick(payload as AnyObj, [
        "trans_prob_final", "transition", "transitions", "M_prob_final", "P", "T",
      ]);
      setTrans(toTransMap(safeParseJSON(transCandidate)));

      const rKey = payload.root_name ?? (pick(payload as AnyObj, ["rootKey", "root_key"]) as string | undefined);
      setRootKey(typeof rKey === "string" ? rKey : null);

      const aPi = as4Vec(payload.D_pi ?? (pick(payload as AnyObj, ["alpha_pi", "dirichlet_pi"]) as unknown)) ?? null;
      const aM  = as4Vec(payload.D_M  ?? (pick(payload as AnyObj, ["alpha_M",  "dirichlet_M"])  as unknown))  ?? null;
      setAlphaPi(aPi);
      setAlphaM(aM);
    } catch (e) {
      setErr(e instanceof Error ? e.message : "Failed to load");
    } finally {
      setLoading(false);
    }
  }, [job, rep]);

  useEffect(() => { if (job) void load(); }, [job, rep, load]);

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
    return `mle_${job || "job"}${repPart}`;
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

  return (
    <div className="p-6 max-w-[1400px] mx-auto space-y-4" ref={plotsRef}>
      {/* Header row */}
      <div className="flex flex-wrap items-end gap-3">
        <h1 className="text-2xl font-bold">
          Maximum likelihood estimate of GMM and density plot of Dirichlet parameters
        </h1>

        <div className="text-sm ml-2">
          Job: <span className="font-mono">{job || "(none)"} </span>
        </div>

        {/* Repetition selector */}
        <label className="ml-4 text-sm flex items-center gap-2">
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

        {rep ? <div className="text-xs text-gray-600">Using rep: <code>{rep}</code></div> : null}

        <button
          type="button"
          onClick={() => void load()}
          className="ml-auto px-4 py-2 rounded bg-black text-white disabled:opacity-50"
          disabled={!job || loading}
          title="Reload"
        >
          {loading ? "Loading…" : "Reload"}
        </button>
      </div>

      {/* Downloads toolbar (always visible; buttons disabled if nothing to export) */}
      <div className="flex flex-wrap items-center gap-2 border rounded p-3 bg-white">
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
          No Dirichlet data available for this job{rep ? ` (rep ${rep})` : ""}.
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
            showBetaCI
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
