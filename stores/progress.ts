// stores/progress.ts
"use client";

import { create } from "zustand";

/* ===================== Types & constants ===================== */

export type Method = "parsimony" | "dirichlet" | "hss";
export type LayerPerc = { bottom: number; middle: number; top: number };
export type Edge = [string, string, number];

export const LAYERS = ["bottom", "middle", "top"] as const;
export type LayerName = (typeof LAYERS)[number];

// For the long-hand angle
export type LayerIdx = 0 | 1 | 2;

export type PaletteMap = Record<Method, Record<LayerName, `#${string}`>>;

export const CAKE_COLORS: PaletteMap = {
  dirichlet: { bottom: "#D58B2E", middle: "#B06608", top: "#7A4606" },
  parsimony: { bottom: "#BF6A4A", middle: "#943F22", top: "#5E2413" },
  hss:       { bottom: "#2F4E8F", middle: "#09255D", top: "#06183D" },
};

export const DEFAULT_COLORS = CAKE_COLORS;

/* ===================== Full input parameters ===================== */

export type PatternCounts = { coarse: number; medium: number; fine: number };
export type EmInput = {
  thr: number;
  reps: number;
  maxIter: number;
  D_pi: [number, number, number, number];
  D_M: [number, number, number, number];
  requireAlphaGe1: boolean;
  includeParsimony: boolean;
  includeHSS: boolean;
  patternCounts: PatternCounts;

  // optional convenience
  useDefaultSeqs: boolean;
  dnaFileName: string | null;
};

/* ===================== Helpers ===================== */

const normalizePerc = (p?: Partial<LayerPerc>): LayerPerc => {
  const prev: LayerPerc = { bottom: 100, middle: 70, top: 40 };
  let middle = Math.max(0, Math.min(100, Number(p?.middle ?? prev.middle)));
  let top = Math.max(0, Math.min(100, Number(p?.top ?? prev.top)));
  if (top > middle) [top, middle] = [middle, top];
  return { bottom: 100, middle, top };
};

const isInternal = (name: string) => /^h_/.test(name);

const countInternalNodesFromEdges = (edges: Edge[]): number => {
  const s = new Set<string>();
  for (const [u, v] of edges) {
    if (isInternal(u)) s.add(u);
    if (isInternal(v)) s.add(v);
  }
  return s.size;
};

const pushUnique = <T,>(arr: T[], v: T) => {
  if (!arr.includes(v)) arr.push(v);
  return arr;
};

/* ===================== Worker message types ===================== */

type ProgressKind =
  | "setup:reset"
  | "setup:reps"
  | "setup:edges"
  | "rep"
  | "method:start"
  | "node"
  | "root:done"
  | "layer:done";

type MsgEmtrifle = {
  type: "emtrifle";
  payload: { method: Method; root: string; rep: number };
};

type MsgEdges = { type: "edges"; edges: Edge[] };

type MsgProgress = {
  type: "progress";
  kind: ProgressKind;
  totalReps?: number;
  nodeTotal?: number;
  rep?: number;
  method?: Method;
  nodeIndex?: number;
  rootName?: string;
  layer?: LayerIdx;
};

type MsgLog = { type: "log"; line: string };

type MsgNewick = { type: "tree_newick"; payload: string };

type WorkerMsg = MsgEmtrifle | MsgEdges | MsgProgress | MsgNewick | MsgLog;

function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null;
}

function isEmtrifleMsg(x: unknown): x is MsgEmtrifle {
  if (!isObject(x) || (x as { type?: unknown }).type !== "emtrifle") return false;
  const p = (x as { payload?: unknown }).payload as Record<string, unknown> | undefined;
  return !!p &&
    (p.method === "parsimony" || p.method === "dirichlet" || p.method === "hss") &&
    typeof p.root === "string" &&
    Number.isFinite((p as { rep?: unknown }).rep as number);
}

function isEdgesMsg(x: unknown): x is MsgEdges {
  return isObject(x) && (x as { type?: unknown }).type === "edges" && Array.isArray((x as { edges?: unknown }).edges);
}

function isProgressMsg(x: unknown): x is MsgProgress {
  if (!isObject(x) || (x as { type?: unknown }).type !== "progress") return false;
  const kind = (x as { kind?: unknown }).kind;
  return (
    kind === "setup:reset" ||
    kind === "setup:reps" ||
    kind === "setup:edges" ||
    kind === "rep" ||
    kind === "method:start" ||
    kind === "node" ||
    kind === "root:done" ||
    kind === "layer:done"
  );
}

function isLogMsg(x: unknown): x is MsgLog {
  return isObject(x) && (x as { type?: unknown }).type === "log" && typeof (x as { line?: unknown }).line === "string";
}

function isNewickMsg(x: unknown): x is MsgNewick {
  return isObject(x) && (x as { type?: unknown }).type === "tree_newick" && typeof (x as { payload?: unknown }).payload === "string";
}

/* ===================== Run inputs (fractions) ===================== */

type RunInputs = {
  reps: number;
  layerPercGlobal?: Partial<LayerPerc>;
  layerPercByMethod?: Partial<Record<Method, Partial<LayerPerc>>>;
};

/* ===================== Store shape ===================== */

type State = {
  // progress / clock (header)
  currentMethod: Method | null;
  currentRoot: string | null;
  nodeIndex: number;
  nodeTotal: number;

  // short hand
  currentRep: number;
  completedReps: number;
  totalReps: number;

  // last completed layer in the *current header*
  lastLayer: LayerIdx | null;

  // sticky long-hand (advance only at end of layer 3)
  longHandMethod: Method | null;
  longHandNodeIndex: number | null;   // 1-based
  longHandNodeTotal: number | null;

  nodeOrder: Record<Method, string[]>;
  completedByNode: Record<string, Method | true>;
  completedNodes: Set<string>;

  // fractions input (for wedge %s)
  inputReps: number | null;
  inputLayerPercGlobal: LayerPerc | null;
  inputLayerPercByMethod: Partial<Record<Method, LayerPerc>>;

  // data sources
  edges: Edge[];
  newick: string | null;

  // palette
  colors: PaletteMap;
  clockSingleColor: `#${string}`;

  // full input parameters (ECDLL)
  emInput: EmInput;

  // logs (for DebugProgressPane)
  logs: string[];

  // internals
  _seenReps: Record<Method, Set<number>>;
  _workerAttached: boolean;
  _worker?: Worker | null;
  _workerOnMsg?: (evt: MessageEvent) => void;
};

type Actions = {
  // fractions-based setup (existing)
  ingestRunInputs: (args: RunInputs) => void;
  setTotalReps: (n: number) => void;
  setNewick: (s: string) => void;

  // quick setters
  setLastLayer: (l: LayerIdx | null) => void;
  setCurrentMethod: (m: Method | null) => void;

  // worker linkage
  attachWorker: (w: Worker | null) => void;
  handleEmtrifle: (payload: { method: Method; root: string; rep: number }) => void;

  // logs
  pushLog: (line: string) => void;
  clearLogs: () => void;

  // readers
  fractionsFor: (m: Method, fallback?: [number, number, number]) => [number, number, number];
  effectiveReps: () => number;

  // palette ops
  setColors: (patch: Partial<Record<Method, Partial<Record<LayerName, `#${string}`>>>>) => void;
  cakeColorFor: (m: Method, layer: LayerName | 0 | 1 | 2) => `#${string}`;
  cakeMethodColors: (m: Method) => [bottom: `#${string}`, middle: `#${string}`, top: `#${string}`];
  setClockSingleColor: (hex: `#${string}`) => void;

  // input params ops
  setEmInput: (patch: Partial<EmInput>) => void;

  // parse raw “[EM_Method] {json}”
  handleProgressLine: (line: string) => void;
};

/* ===================== Store ===================== */

export const useProgressStore = create<State & Actions>((set, get) => ({
  // progress defaults
  currentMethod: null,
  currentRoot: null,
  nodeIndex: 0,
  nodeTotal: 0,

  currentRep: 0,
  completedReps: 0,
  totalReps: 30, // sensible default

  // header layer marker (not used for long-hand anymore)
  lastLayer: null,

  // sticky long-hand defaults
  longHandMethod: null,
  longHandNodeIndex: null,
  longHandNodeTotal: null,

  nodeOrder: { parsimony: [], dirichlet: [], hss: [] },

  completedByNode: {},
  completedNodes: new Set<string>(),

  inputReps: null,
  inputLayerPercGlobal: null,
  inputLayerPercByMethod: {},

  edges: [],
  newick: null,

  colors: CAKE_COLORS,
  clockSingleColor: "#111111",

  // default input-form parameters
  emInput: {
    thr: 0.01,
    reps: 30,
    maxIter: 100,
    D_pi: [100, 100, 100, 100],
    D_M: [100, 2, 2, 2],
    requireAlphaGe1: true,
    includeParsimony: true,
    includeHSS: true,
    patternCounts: { coarse: 26, medium: 42, fine: 100 },
    useDefaultSeqs: false,
    dnaFileName: null,
  },

  // logs
  logs: [],

  _seenReps: {
    parsimony: new Set<number>(),
    dirichlet: new Set<number>(),
    hss: new Set<number>(),
  },
  _workerAttached: false,
  _worker: null,
  _workerOnMsg: undefined,

  /* ---------------- Actions ---------------- */

  ingestRunInputs: (args) => {
    const reps = Math.max(1, Number(args.reps) || 1);
    const global = args.layerPercGlobal ? normalizePerc(args.layerPercGlobal) : null;

    const perMethod: Partial<Record<Method, LayerPerc>> = {};
    if (args.layerPercByMethod) {
      (["parsimony", "dirichlet", "hss"] as const).forEach((k) => {
        const p = args.layerPercByMethod?.[k];
        if (p) perMethod[k] = normalizePerc(p);
      });
    }

    set({
      inputReps: reps,
      totalReps: reps,
      inputLayerPercGlobal: global,
      inputLayerPercByMethod: perMethod,
    });
  },

  setTotalReps: (n) => set({ totalReps: Math.max(1, Number(n) || 1) }),
  setNewick: (s) => set({ newick: s || null }),

  setLastLayer: (l) => set({ lastLayer: l }),
  setCurrentMethod: (m) => set({ currentMethod: m }),

  handleEmtrifle: ({ method, root, rep }) => {
    if (!method || !root) return;

    const { nodeOrder, completedNodes, _seenReps, completedByNode } = get();

    const arr = nodeOrder[method].slice();
    pushUnique(arr, root);

    const nextDone = new Set(completedNodes);
    nextDone.add(root);
    const nextCompletedByNode = { ...completedByNode, [root]: method };

    _seenReps[method].add(rep);
    const globalReps = new Set<number>([
      ..._seenReps.parsimony,
      ..._seenReps.dirichlet,
      ..._seenReps.hss,
    ]).size;

    set({
      currentMethod: method,
      currentRoot: root,
      nodeOrder: { ...nodeOrder, [method]: arr },
      completedNodes: nextDone,
      completedByNode: nextCompletedByNode,
      completedReps: globalReps,
      // NOTE: do NOT set currentRep here; only "rep" messages should move the short hand
    });
  },

  // Logs: keep a rolling buffer
  pushLog: (line) =>
    set((s) => {
      const MAX = 2000; // keep last 2k lines
      const trimmed = String(line ?? "").trim();
      if (!trimmed) return {};
      const next = s.logs.length >= MAX ? s.logs.slice(-Math.floor(MAX * 0.9)) : s.logs.slice();
      next.push(trimmed);
      return { logs: next };
    }),

  clearLogs: () => set({ logs: [] }),

  attachWorker: (w) => {
    if (!w) return;

    // Detach any previous worker/listener
    const prev = get()._worker;
    const prevHandler = get()._workerOnMsg;
    if (prev && prevHandler) {
      try { prev.removeEventListener("message", prevHandler); } catch {}
    }

    const onMsg = (evt: MessageEvent) => {
      const msg: unknown = evt.data;

      // Log messages: tee to store + parse progress if applicable
      if (isLogMsg(msg)) {
        const line = (msg.line || "").trim();
        get().pushLog(line);
        get().handleProgressLine(line);
        return;
      }

      if (isEmtrifleMsg(msg)) {
        const { method, root, rep } = msg.payload;
        get().handleEmtrifle({ method, root, rep });
        return;
      }

      if (isEdgesMsg(msg)) {
        const edges = msg.edges;
        const nodeTotal = countInternalNodesFromEdges(edges);
        set({ edges, nodeTotal: nodeTotal || get().nodeTotal });
        return;
      }

      if (isProgressMsg(msg)) {
        switch (msg.kind) {
          case "setup:reset": {
            set({
              currentMethod: null,
              currentRoot: null,
              nodeIndex: 0,
              nodeTotal: 0,
              currentRep: 0,
              completedReps: 0,
              lastLayer: null,
              longHandMethod: null,
              longHandNodeIndex: null,
              longHandNodeTotal: null,
            });
            break;
          }

          case "setup:reps":
            if (Number.isFinite(msg.totalReps)) {
              set({ totalReps: Number(msg.totalReps), lastLayer: null });
            }
            break;

          case "setup:edges":
            if (Number.isFinite(msg.nodeTotal)) set({ nodeTotal: Number(msg.nodeTotal) });
            break;

          case "rep":
            // ONLY here should currentRep advance (end of repetition)
            if (Number.isFinite(msg.rep)) set({ currentRep: Number(msg.rep) });
            break;

          case "method:start":
            // reset header layer marker so the pane knows a new method started
            if (msg.method) set({ currentMethod: msg.method, lastLayer: null });
            break;

          case "node": {
            // update header; do not touch sticky long hand
            const patch: Partial<State> = { lastLayer: null };
            if (Number.isFinite(msg.nodeIndex)) patch.nodeIndex = Number(msg.nodeIndex);
            if (Number.isFinite(msg.nodeTotal)) patch.nodeTotal = Number(msg.nodeTotal);
            if (msg.rootName) patch.currentRoot = msg.rootName;
            set(patch);
            break;
          }

          case "layer:done": {
            // Always keep “header” fresh…
            const header: Partial<State> = {};
            if (msg.method) header.currentMethod = msg.method;
            if (Number.isFinite(msg.nodeIndex)) header.nodeIndex = Number(msg.nodeIndex);
            if (Number.isFinite(msg.nodeTotal)) header.nodeTotal = Number(msg.nodeTotal);
            if (msg.rootName) header.currentRoot = msg.rootName;

            // …but ONLY move the long hand when layer 3 (index 2) completes.
            if (msg.layer === 2) {
              set({
                ...header,
                lastLayer: 2,
                longHandMethod: msg.method ?? null,
                longHandNodeIndex: Number.isFinite(msg.nodeIndex)
                  ? Number(msg.nodeIndex)
                  : get().longHandNodeIndex,
                longHandNodeTotal: Number.isFinite(msg.nodeTotal)
                  ? Number(msg.nodeTotal)
                  : get().longHandNodeTotal,
              });
            } else {
              // Layers 0/1: update header only; keep long hand where it is
              set(header);
            }
            // IMPORTANT: do NOT set currentRep here
            break;
          }

          case "root:done":
            // no-op (reserved)
            break;
        }
        return;
      }

      if (isNewickMsg(msg)) {
        set({ newick: msg.payload });
        return;
      }
    };

    w.addEventListener("message", onMsg);
    set({ _workerAttached: true, _worker: w, _workerOnMsg: onMsg });
  },

  // readers
  fractionsFor: (m, fallback = [1, 0.7, 0.4]) => {
    const s = get();
    const src = s.inputLayerPercByMethod[m] ?? s.inputLayerPercGlobal;
    if (!src) return fallback as [number, number, number];
    return [1, src.middle / 100, src.top / 100];
  },

  effectiveReps: () => {
    const s = get();
    return s.inputReps ?? s.totalReps ?? 1;
  },

  // palette ops
  setColors: (patch) =>
    set((s) => {
      const next = { ...s.colors };
      (Object.keys(patch) as Method[]).forEach((m) => {
        next[m] = { ...next[m], ...(patch[m] || {}) };
      });
      return { colors: next };
    }),

  cakeColorFor: (m, layer) => {
    const s = get();
    const key = typeof layer === "number" ? LAYERS[Math.max(0, Math.min(2, layer))] : layer;
    return s.colors[m][key];
  },

  cakeMethodColors: (m) => {
    const s = get();
    const c = s.colors[m];
    return [c.bottom, c.middle, c.top];
  },

  setClockSingleColor: (hex) => set({ clockSingleColor: hex }),

  // store/update all input parameters; if reps present, also sync totalReps
  setEmInput: (patch) =>
    set((s) => {
      const emInput: EmInput = { ...s.emInput, ...patch };
      const next: Partial<State> = { emInput };
      if (typeof patch.reps === "number" && Number.isFinite(patch.reps)) {
        next.totalReps = Math.max(1, patch.reps);
      }
      return next as State;
    }),

  // Accept lines like: "[EM_Method] {json}"
  handleProgressLine: (line: string) => {
    const trimmed = (line || "").trim();
    const i = trimmed.indexOf("{");
    if (i < 0) return;
    try {
      const j = JSON.parse(trimmed.slice(i));
      if (j && j.type === "progress") {
        const msg: MsgProgress = { type: "progress", ...(j as Record<string, unknown>) } as MsgProgress;
        if (!isProgressMsg(msg)) return;

        switch (msg.kind) {
          case "setup:reset":
            set({
              currentMethod: null,
              currentRoot: null,
              nodeIndex: 0,
              nodeTotal: 0,
              currentRep: 0,
              completedReps: 0,
              lastLayer: null,
              longHandMethod: null,
              longHandNodeIndex: null,
              longHandNodeTotal: null,
            });
            break;

          case "setup:reps":
            if (Number.isFinite(msg.totalReps)) set({ totalReps: Number(msg.totalReps), lastLayer: null });
            break;

          case "setup:edges":
            if (Number.isFinite(msg.nodeTotal)) set({ nodeTotal: Number(msg.nodeTotal) });
            break;

          case "rep":
            // ONLY here should currentRep advance (end of repetition)
            if (Number.isFinite(msg.rep)) set({ currentRep: Number(msg.rep) });
            break;

          case "method:start":
            if (msg.method) set({ currentMethod: msg.method, lastLayer: null });
            break;

          case "node": {
            const patch: Partial<State> = { lastLayer: null };
            if (Number.isFinite(msg.nodeIndex)) patch.nodeIndex = Number(msg.nodeIndex);
            if (Number.isFinite(msg.nodeTotal)) patch.nodeTotal = Number(msg.nodeTotal);
            if (msg.rootName) patch.currentRoot = msg.rootName;
            set(patch);
            break;
          }

          case "layer:done": {
            // Always keep “header” fresh…
            const header: Partial<State> = {};
            if (msg.method) header.currentMethod = msg.method;
            if (Number.isFinite(msg.nodeIndex)) header.nodeIndex = Number(msg.nodeIndex);
            if (Number.isFinite(msg.nodeTotal)) header.nodeTotal = Number(msg.nodeTotal);
            if (msg.rootName) header.currentRoot = msg.rootName;

            // …but ONLY move the long hand when layer 3 (index 2) completes.
            if (msg.layer === 2) {
              set({
                ...header,
                lastLayer: 2,
                longHandMethod: msg.method ?? null,
                longHandNodeIndex: Number.isFinite(msg.nodeIndex)
                  ? Number(msg.nodeIndex)
                  : get().longHandNodeIndex,
                longHandNodeTotal: Number.isFinite(msg.nodeTotal)
                  ? Number(msg.nodeTotal)
                  : get().longHandNodeTotal,
              });
            } else {
              // Layers 0/1: update header only; keep long hand where it is
              set(header);
            }
            // IMPORTANT: do NOT set currentRep here
            break;
          }

          case "root:done":
            // reserved
            break;
        }
      }
    } catch {
      // ignore parse errors
    }
  },
}));
