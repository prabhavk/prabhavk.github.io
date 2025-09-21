// components/DaylightTree.tsx
"use client";

import React from "react";
import type { Method } from "@/stores/progress";

/** -------- Edge input (both tuple and object accepted) -------- */
type EdgeInput =
  | [string, string, number]
  | { source: string; target: string; length?: number | null };

type Props = {
  className?: string;
  edges: ReadonlyArray<EdgeInput>;
  /** (ignored by layout; we root at a pre-existing internal node nearest the tree center) */
  currentRoot?: string;
  /** method whose color to use for live/incomplete nodes */
  activeMethod?: Method | null;
  /** per-node completion; true = completed for all methods, or a Method = completed for that method only */
  completedByNode: Record<string, true | Method>;
  /** exactly one internal node that is currently active (running) */
  activeNode?: string | null;
  padding?: number;
  /** Show both leaf and internal-node labels */
  showLabels?: boolean;
  /** Number of global passes for the daylight algorithm (default 1) */
  daylightPasses?: number;
  /** Visual magnification (applied around canvas center). Default 1.5x for readability */
  zoom?: number;
  /** Vertical pixel shift applied outside the zoom transform (moves the whole tree up/down) */
  yShift?: number;
};

/* ===== Colors ===== */
const METHOD_COLOR: Record<Method, string> = {
  parsimony: "#10b981",
  dirichlet: "#f59e0b",
  hss: "#8b5cf6",
};
const NEUTRAL = "#64748b";

const TAU = Math.PI * 2;
const PI = Math.PI;

/* ---------- helpers ---------- */
function toTuple(e: EdgeInput): [string, string, number] {
  if (Array.isArray(e)) {
    const [u, v, w] = e;
    return [u, v, Number.isFinite(w) ? Number(w) : 1];
  }
  return [e.source, e.target, Number.isFinite(e.length ?? NaN) ? Number(e.length) : 1];
}
function wrap01(x: number) {
  const r = x - Math.floor(x);
  return r < 0 ? r + 1 : r;
}
function as01Angle(a: number) {
  return wrap01(a / TAU) * TAU;
}
function polar(cx: number, cy: number, r: number, theta: number) {
  // 12 o'clock = 0
  const x = cx + Math.cos(theta - PI / 2) * r;
  const y = cy + Math.sin(theta - PI / 2) * r;
  return { x, y };
}

/* ---------- graph ---------- */
type Adj = Map<string, Array<{ to: string; w: number }>>;
function buildAdj(edges: ReadonlyArray<EdgeInput>): Adj {
  const adj: Adj = new Map();
  for (const e of edges) {
    const [u, v, w] = toTuple(e);
    if (!adj.has(u)) adj.set(u, []);
    if (!adj.has(v)) adj.set(v, []);
    adj.get(u)!.push({ to: v, w });
    adj.get(v)!.push({ to: u, w });
  }
  return adj;
}
function degrees(adj: Adj): Map<string, number> {
  const deg = new Map<string, number>();
  for (const [u, list] of adj.entries()) {
    const set = new Set<string>();
    for (const { to } of list) set.add(to);
    deg.set(u, set.size);
  }
  return deg;
}
function getEdgeWeight(adj: Adj, u: string, v: string): number {
  const hit = adj.get(u)?.find((e) => e.to === v);
  return hit ? hit.w : 1;
}

/* ---------- diameter midpoint -> choose real node root ---------- */
type FarthestResult = {
  dist: Map<string, number>;
  parent: Map<string, string | null>;
  farNode: string;
};
function farthestFrom(adj: Adj, start: string): FarthestResult {
  const dist = new Map<string, number>();
  const parent = new Map<string, string | null>();
  const stack: Array<{ u: string; p: string | null; d: number }> = [{ u: start, p: null, d: 0 }];
  dist.set(start, 0);
  parent.set(start, null);
  while (stack.length) {
    const { u, p, d } = stack.pop()!;
    for (const { to: v, w } of adj.get(u) || []) {
      if (v === p) continue;
      dist.set(v, d + w);
      parent.set(v, u);
      stack.push({ u: v, p: u, d: d + w });
    }
  }
  let farNode = start;
  let maxD = -Infinity;
  for (const [k, di] of dist.entries()) {
    if (di > maxD) { maxD = di; farNode = k; }
  }
  return { dist, parent, farNode };
}
function pathBetween(parent: Map<string, string | null>, src: string, dst: string) {
  const path: string[] = [];
  let cur: string | null = dst;
  while (cur !== null) {
    path.push(cur);
    if (cur === src) break;
    cur = parent.get(cur) ?? null;
  }
  path.reverse();
  return path;
}
function pickCentralRoot(adj: Adj, deg: Map<string, number>): string {
  const first = adj.keys().next();
  if (first.done) return "";
  const a = first.value;
  const A = farthestFrom(adj, a);
  const B = farthestFrom(adj, A.farNode);
  const src = B.farNode;
  const C = farthestFrom(adj, src);
  const dst = C.farNode;

  const { parent } = farthestFrom(adj, src);
  const path = pathBetween(parent, src, dst);

  const cum: number[] = [0];
  for (let i = 1; i < path.length; i++) {
    const u = path[i - 1], v = path[i];
    const w = (adj.get(u)!.find((x) => x.to === v)!.w);
    cum.push(cum[i - 1] + w);
  }
  const total = cum[cum.length - 1];
  const mid = total / 2;

  let best = 0, bestAbs = Math.abs(cum[0] - mid);
  for (let i = 1; i < cum.length; i++) {
    const d = Math.abs(cum[i] - mid);
    if (d < bestAbs) { bestAbs = d; best = i; }
  }
  if ((deg.get(path[best]) || 0) <= 1) {
    let L = best - 1, R = best + 1;
    while (L >= 0 || R < path.length) {
      if (L >= 0 && (deg.get(path[L]) || 0) > 1) return path[L];
      if (R < path.length && (deg.get(path[R]) || 0) > 1) return path[R];
      L--; R++;
    }
  }
  return path[best];
}

/* ---------- root at chosen node ---------- */
type Rooted = {
  root: string;
  nodes: string[];
  parent: Map<string, string | null>;
  children: Map<string, string[]>;
  deg: Map<string, number>;
  adj: Adj;
};
function rootAtNode(adj: Adj, rootId: string): Rooted {
  if (!rootId) return { root: "", nodes: [], parent: new Map(), children: new Map(), deg: new Map(), adj };

  const parent = new Map<string, string | null>();
  const children = new Map<string, string[]>();
  const nodes = Array.from(adj.keys());
  for (const u of nodes) children.set(u, []);

  const stack: Array<{ u: string; p: string | null }> = [{ u: rootId, p: null }];
  parent.set(rootId, null);
  while (stack.length) {
    const { u, p } = stack.pop()!;
    for (const { to: v } of adj.get(u) || []) {
      if (v === p) continue;
      parent.set(v, u);
      children.get(u)!.push(v);
      stack.push({ u: v, p: u });
    }
  }
  const deg = degrees(adj);
  return { root: rootId, nodes: Array.from(children.keys()), parent, children, deg, adj };
}

/* ---------- equal-angle layout ---------- */
type LayoutEA = {
  angle: Map<string, number>;
  depth: Map<string, number>;
  maxDepth: number;
  leaves: string[];
};
function buildEqualAngleLayout(
  root: string,
  children: Map<string, string[]>,
  deg: Map<string, number>
): LayoutEA {
  const leaves: string[] = [];
  (function dfsLeaves(u: string) {
    const kids = children.get(u) || [];
    if (!kids.length) { leaves.push(u); return; }
    for (const v of kids) dfsLeaves(v);
  })(root);
  const L = Math.max(1, leaves.length);

  const angle = new Map<string, number>();
  leaves.forEach((u, i) => {
    const a = TAU * (i / L);
    angle.set(u, as01Angle(a));
  });

  const postorder: string[] = [];
  (function buildPost(u: string) {
    for (const v of children.get(u) || []) buildPost(v);
    postorder.push(u);
  })(root);

  const leafCount = new Map<string, number>();
  for (const u of postorder) {
    const kids = children.get(u) || [];
    if (!kids.length) { leafCount.set(u, 1); continue; }
    let cnt = 0, sx = 0, sy = 0;
    for (const v of kids) {
      const c = leafCount.get(v) || 0;
      cnt += c;
      const av = angle.get(v);
      if (av != null) { sx += c * Math.cos(av); sy += c * Math.sin(av); }
    }
    leafCount.set(u, cnt);
    if (cnt > 0) {
      const a = Math.atan2(sy, sx);
      angle.set(u, as01Angle(a));
    }
  }
  if (angle.get(root) == null) angle.set(root, 0);

  const depth = new Map<string, number>();
  let maxDepth = 0;
  (function dfsDepth(u: string, d: number) {
    depth.set(u, d);
    maxDepth = Math.max(maxDepth, d);
    for (const v of children.get(u) || []) dfsDepth(v, d + 1);
  })(root, 0);

  return { angle, depth, maxDepth: Math.max(1, maxDepth), leaves };
}

/* ---------- Euler tour helpers ---------- */
type EulerInfo = {
  tin: Map<string, number>;
  tout: Map<string, number>;
};
function eulerize(root: string, children: Map<string, string[]>): EulerInfo {
  const tin = new Map<string, number>();
  const tout = new Map<string, number>();
  let T = 0;
  const stack: Array<{ u: string; i: number }> = [{ u: root, i: 0 }];
  while (stack.length) {
    const top = stack[stack.length - 1];
    if (top.i === 0) tin.set(top.u, T++);
    const kids = children.get(top.u) || [];
    if (top.i < kids.length) {
      const v = kids[top.i++];
      stack.push({ u: v, i: 0 });
    } else {
      tout.set(top.u, T);
      stack.pop();
    }
  }
  return { tin, tout };
}
function isInSubtree(tin: Map<string, number>, tout: Map<string, number>, a: string, b: string) {
  const ta = tin.get(a)!;
  const tb = tin.get(b)!;
  const te = tout.get(b)!;
  return tb <= ta && ta < te;
}

/* =================== preorder =================== */
export function buildPreorderInternalList(
  root: string,
  children: Map<string, string[]>,
  deg: Map<string, number>
): string[] {
  const out: string[] = [];
  if (!root) return out;
  const stack: string[] = [root];
  while (stack.length) {
    const u = stack.pop()!;
    if ((deg.get(u) || 0) > 1) out.push(u);
    const kids = children.get(u);
    if (kids && kids.length) {
      for (let i = kids.length - 1; i >= 0; i--) stack.push(kids[i]);
    }
  }
  return out;
}

export function computeSubtreeWedgesAtNode(
  p: string,
  params: {
    parent: Map<string, string | null>;
    children: Map<string, string[]>;
    leaves: string[];
    angle: Map<string, number>;
    depth: Map<string, number>;
    radiusOfDepth: (d: number) => number;
    cx: number;
    cy: number;
    tin: Map<string, number>;
    tout: Map<string, number>;
  }
) {
  const { parent, children, leaves, angle, depth, radiusOfDepth, cx, cy, tin, tout } = params;

  const nodePos = (u: string) => {
    const a = angle.get(u) ?? 0;
    const d = depth.get(u) ?? 0;
    return polar(cx, cy, radiusOfDepth(d), a);
  };

  const rayAngleFromP = (q: string) => {
    const P = nodePos(p);
    const Q = nodePos(q);
    return as01Angle(Math.atan2(Q.y - P.y, Q.x - P.x));
  };

  const nbrs: Array<{ kind: "child" | "parent"; id?: string; dir: number }> = [];
  for (const c of (children.get(p) || [])) {
    nbrs.push({ kind: "child", id: c, dir: rayAngleFromP(c) });
  }
  const par = parent.get(p);
  if (par !== null && par !== undefined) {
    nbrs.push({ kind: "parent", id: par, dir: rayAngleFromP(par) });
  }

  nbrs.sort((a, b) => a.dir - b.dir);

  const wedges: Array<{
    kind: "child" | "parent";
    id?: string;
    left: number;
    right: number;
    center: number;
    width: number;
  }> = [];
  const P = nodePos(p);

  const raysToLeaves = (leafIds: string[]) => {
    const arr: number[] = [];
    for (const L of leafIds) {
      const Lp = nodePos(L);
      arr.push(as01Angle(Math.atan2(Lp.y - P.y, Lp.x - P.x)));
    }
    return arr;
  };

  for (const nb of nbrs) {
    let leafAngles: number[] = [];
    if (nb.kind === "child" && nb.id) {
      const ls = leaves.filter((L) => isInSubtree(tin, tout, L, nb.id!));
      leafAngles = raysToLeaves(ls);
    } else {
      const ls = leaves.filter((L) => !isInSubtree(tin, tout, L, p));
      leafAngles = raysToLeaves(ls);
    }

    if (!leafAngles.length) leafAngles = [nb.dir];

    const shifted = leafAngles
      .map((a) => {
        let d = a - nb.dir;
        while (d <= -PI) d += TAU;
        while (d > PI) d -= TAU;
        return nb.dir + d;
      })
      .sort((x, y) => x - y);

    const left = shifted[0];
    const right = shifted[shifted.length - 1];
    const width = Math.max(0, right - left);
    const center = (left + right) / 2;

    wedges.push({
      kind: nb.kind,
      id: nb.id,
      left: as01Angle(left),
      right: as01Angle(right),
      center: as01Angle(center),
      width,
    });
  }

  return wedges;
}

/* ==================== daylight helpers & rotations ==================== */

export function applyDeltaToSubtree(
  subRoot: string,
  delta: number,
  angle: Map<string, number>,
  tin: Map<string, number>,
  tout: Map<string, number>
) {
  if (!Number.isFinite(delta) || delta === 0) return;
  const t0 = tin.get(subRoot)!;
  const t1 = tout.get(subRoot)!;
  for (const [id, a] of angle.entries()) {
    const ti = tin.get(id);
    if (ti != null && ti >= t0 && ti < t1) {
      let v = a + delta;
      v = ((v % TAU) + TAU) % TAU;
      angle.set(id, v);
    }
  }
}

function unwrapTriple(
  a: { left: number; right: number },
  b: { left: number; right: number },
  c: { left: number; right: number }
) {
  const rA = 0;
  const lA = ((a.left - a.right + TAU) + TAU) % TAU;

  const lB = ((b.left - a.right + TAU) + TAU) % TAU;
  let rB = ((b.right - a.right + TAU) + TAU) % TAU;
  let lC = ((c.left - a.right + TAU) + TAU) % TAU;
  let rC = ((c.right - a.right + TAU) + TAU) % TAU;

  if (rB < lB) rB += TAU;
  if (lC < rB) { lC += TAU; rC += TAU; }
  if (rC < lC) rC += TAU;

  return { rA, lA, lB, rB, lC, rC };
}

export function equalizeDaylightAtNode_Generic(
  u: string,
  visitedInternals: Set<string>,
  params: {
    parent: Map<string, string | null>;
    children: Map<string, string[]>;
    leaves: string[];
    angle: Map<string, number>;
    depth: Map<string, number>;
    radiusOfDepth: (d: number) => number;
    cx: number;
    cy: number;
    tin: Map<string, number>;
    tout: Map<string, number>;
  }
): { deltaB: number; deltaC: number } | null {
  const { children, angle, tin, tout } = params;

  const wedges = computeSubtreeWedgesAtNode(u, params);
  if (wedges.length < 3) return null;

  const three =
  wedges.length === 3 ? wedges : [...wedges].sort((x, y) => y.width - x.width).slice(0, 3);

  const isLeafChild = (w: { kind: "child" | "parent"; id?: string }) =>
    w.kind === "child" && w.id ? (children.get(w.id!)?.length ?? 0) === 0 : false;

  const visitedIdxs = three
    .map((w, i) => (w.id && visitedInternals.has(w.id) ? i : -1))
    .filter((i) => i >= 0);

  let fixedIdx: number;
  if (visitedIdxs.length) {
    fixedIdx = visitedIdxs[0]!;
  } else {
    const nonLeafIdxs = [0, 1, 2].filter((i) => !isLeafChild(three[i]!));
    fixedIdx = nonLeafIdxs.length
      ? nonLeafIdxs.reduce(
          (best, i) => (three[i]!.width > three[best]!.width ? i : best),
          nonLeafIdxs[0]!
        )
      : 0;
  }

  const a = three[fixedIdx]!;
  const b = three[(fixedIdx + 1) % 3]!;
  const c = three[(fixedIdx + 2) % 3]!;

  console.debug(
    `[Daylight] Node=${u} fixed=${a.id ?? a.kind} rotated=[${b.id ?? b.kind}, ${c.id ?? c.kind}]`
  );

  const Wa = ((a.right - a.left + TAU) + TAU) % TAU;
  const Wb = ((b.right - b.left + TAU) + TAU) % TAU;
  const Wc = ((c.right - c.left + TAU) + TAU) % TAU;
  const G = (TAU - (Wa + Wb + Wc)) / 3;

  const { rA, lA, lB, rB, lC, rC } = unwrapTriple(a, b, c);
  const gap_ab = lB - rA;
  const gap_ca = (lA + TAU) - rC;

  const deltaB = G - gap_ab;
  const deltaC = gap_ca - G;

  if (b.id) applyDeltaToSubtree(b.id, deltaB, angle, tin, tout);
  if (c.id) applyDeltaToSubtree(c.id, deltaC, angle, tin, tout);

  return { deltaB, deltaC };
}

function angleDistance(a: number, b: number): number {
  const d = Math.abs(as01Angle(a) - as01Angle(b));
  return Math.min(d, TAU - d);
}

export function daylightAngleBetween(A: { left: number; right: number }, B: { left: number; right: number }) {
  const d1 = angleDistance(A.left, B.right);
  const d2 = angleDistance(A.right, B.left);
  return Math.min(d1, d2);
}

/* ================================================================
   =====================  COMPONENT (Daylight)  ====================
   ================================================================ */

export default function DaylightTree({
  className = "",
  edges,
  activeMethod = null,
  completedByNode,
  activeNode = null,
  padding = 12,
  showLabels = true,
  daylightPasses = 1,
  zoom = 0.9,
  yShift = -40,
}: Props) {
  const ref = React.useRef<HTMLDivElement>(null);
  const [size, setSize] = React.useState({ w: 0, h: 0 });

  React.useEffect(() => {
    if (!ref.current || typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver((entries) => {
      for (const e of entries) setSize({ w: e.contentRect.width, h: e.contentRect.height });
    });
    ro.observe(ref.current);
    return () => ro.disconnect();
  }, []);

  /* ---------- Build static structure from edges ---------- */
  const core = React.useMemo(() => {
    if (!edges?.length) {
      return {
        nodes: [] as string[],
        root: "",
        parent: new Map<string, string | null>(),
        children: new Map<string, string[]>(),
        deg: new Map<string, number>(),
        adj: new Map() as Adj,
      };
    }
    const adj = buildAdj(edges);
    const deg = degrees(adj);
    const rootNode = pickCentralRoot(adj, deg);
    const rooted = rootAtNode(adj, rootNode);
    return {
      nodes: rooted.nodes,
      root: rooted.root,
      parent: rooted.parent,
      children: rooted.children,
      deg: rooted.deg,
      adj: rooted.adj,
    };
  }, [edges]);

  /* ---------- Equal-angle baseline ---------- */
  const equalAngle = React.useMemo(() => {
    if (!core.root) {
      return {
        angle: new Map<string, number>(),
        depth: new Map<string, number>(),
        maxDepth: 1,
        leaves: [] as string[],
      };
    }
    return buildEqualAngleLayout(core.root, core.children, core.deg);
  }, [core.root, core.children, core.deg]);

  /* ---------- Weighted radii (preserve branch lengths) ---------- */
  const { weightedDepth, weightedMaxDepth } = React.useMemo(() => {
    const dist = new Map<string, number>();
    if (!core.root) return { weightedDepth: dist, weightedMaxDepth: 1 };

    const stack: string[] = [core.root];
    dist.set(core.root, 0);
    while (stack.length) {
      const u = stack.pop()!;
      const base = dist.get(u)!;
      for (const v of core.children.get(u) || []) {
        const w = getEdgeWeight(core.adj, u, v);
        dist.set(v, base + (Number.isFinite(w) ? w : 1));
        stack.push(v);
      }
    }
    let max = 0;
    for (const d of dist.values()) max = Math.max(max, d);
    return { weightedDepth: dist, weightedMaxDepth: Math.max(1e-6, max) };
  }, [core.root, core.children, core.adj]);

  /* ---------- Radii / canvas ---------- */
  const w = Math.max(0, size.w);
  const h = Math.max(0, size.h);
  const cx = w / 2;
  const cy = h / 2;

  const { ringOuter, ringInner } = React.useMemo(() => {
    const R = Math.max(0, Math.min(w, h) / 2 - padding);
    return { ringOuter: R, ringInner: Math.max(0, R * 0.12) };
  }, [w, h, padding]);

  const radiusOfDepth = React.useCallback(
    (d: number) =>
      d === 0 ? 0 : ringInner + (ringOuter - ringInner) * (d / weightedMaxDepth),
    [ringInner, ringOuter, weightedMaxDepth]
  );

  // Euler tour (for subtree rotations)
  const { tin, tout } = React.useMemo(
    () => (core.root ? eulerize(core.root, core.children) : { tin: new Map(), tout: new Map() }),
    [core.root, core.children]
  );

  // Internal nodes in preorder (pass order)
  const internalPre = React.useMemo(
    () => (core.root ? buildPreorderInternalList(core.root, core.children, core.deg) : []),
    [core.root, core.children, core.deg]
  );

  /* ---------- Daylight-adjusted angles + visited internals (N PASSES) ---------- */
  const daylight = React.useMemo(() => {
    const angle = new Map(equalAngle.angle);
    const visitedForColor = new Set<string>();

    if (!core.root || internalPre.length === 0) {
      return { angle, visitedInternals: visitedForColor };
    }

    const passes = Math.max(1, Math.floor(daylightPasses));
    for (let pass = 0; pass < passes; pass++) {
      const visitedThisPass = new Set<string>();
      for (const u of internalPre) {
        equalizeDaylightAtNode_Generic(u, visitedThisPass, {
          parent: core.parent,
          children: core.children,
          leaves: equalAngle.leaves,
          angle,
          depth: weightedDepth,
          radiusOfDepth,
          cx,
          cy,
          tin,
          tout,
        });
        visitedThisPass.add(u);
        visitedForColor.add(u);
      }
      console.debug(`[Daylight] Completed pass ${pass + 1}/${passes}`);
    }

    return { angle, visitedInternals: visitedForColor };
  }, [
    core.root,
    core.parent,
    core.children,
    core.deg,
    equalAngle.leaves,
    internalPre,
    radiusOfDepth,
    cx,
    cy,
    tin,
    tout,
    equalAngle.angle,
    daylightPasses,
    weightedDepth,
  ]);

  // ===== Daylight imbalance per internal node =====
  const daylightImbalance = React.useMemo(() => {
    if (!core.root) return [] as Array<{ id: string; imbalanceRad: number }>;
    const out: Array<{ id: string; imbalanceRad: number }> = [];

    for (const u of core.nodes) {
      if ((core.deg.get(u) || 0) <= 1) continue; // leaves only
      const wedgesAll = computeSubtreeWedgesAtNode(u, {
        parent: core.parent,
        children: core.children,
        leaves: equalAngle.leaves,
        angle: daylight.angle,
        depth: weightedDepth,
        radiusOfDepth,
        cx,
        cy,
        tin,
        tout,
      });

      if (wedgesAll.length < 3) continue;

      const three =
        wedgesAll.length === 3 ? wedgesAll : [...wedgesAll].sort((a, b) => b.width - a.width).slice(0, 3);

      const ordered = [...three].sort((a, b) => a.center - b.center);
      const gaps = [
        daylightAngleBetween(ordered[0], ordered[1]),
        daylightAngleBetween(ordered[1], ordered[2]),
        daylightAngleBetween(ordered[2], ordered[0]),
      ];
      const maxGap = Math.max(...gaps);
      const minGap = Math.min(...gaps);
      out.push({ id: u, imbalanceRad: maxGap - minGap });
    }

    out.sort((a, b) => b.imbalanceRad - a.imbalanceRad);
    return out;
  }, [
    core.root,
    core.nodes,
    core.deg,
    core.parent,
    core.children,
    equalAngle.leaves,
    daylight.angle,
    weightedDepth,
    radiusOfDepth,
    cx,
    cy,
    tin,
    tout,
  ]);

  // Live color always follows the active method (or neutral).
  const liveColor = React.useMemo(
    () => (activeMethod ? METHOD_COLOR[activeMethod] : NEUTRAL),
    [activeMethod]
  );

  // Links — now rendered from DAYLIGHT angles (with weighted radii)
  const links = React.useMemo(() => {
    const out: Array<{ x1: number; y1: number; x2: number; y2: number; key: string }> = [];
    for (const [v, p] of core.parent.entries()) {
      if (p === null) continue;
      const av = daylight.angle.get(v) ?? 0;
      const ap = daylight.angle.get(p) ?? 0;

      const dv = (weightedDepth.get(v) ?? 0);
      const dp = (weightedDepth.get(p) ?? 0);

      const pv = polar(cx, cy, radiusOfDepth(dv), av);
      const pp = polar(cx, cy, radiusOfDepth(dp), ap);
      out.push({ x1: pv.x, y1: pv.y, x2: pp.x, y2: pp.y, key: `${p}->${v}` });
    }
    return out;
  }, [core.parent, daylight.angle, cx, cy, radiusOfDepth, weightedDepth]);

  const nodesViz = React.useMemo(() => {
    const out: Array<{
      id: string;
      x: number;
      y: number;
      rOuter: number;
      color: string;
      showWhiteCap: boolean;
      showPlus: boolean;
      isLeaf: boolean;
    }> = [];

    for (const u of core.nodes) {
      const a = daylight.angle.get(u) ?? 0;
      const d = weightedDepth.get(u) ?? 0;
      const { x, y } = polar(cx, cy, radiusOfDepth(d), a);

      const isLeaf = (core.deg.get(u) || 0) <= 1;
      const isInternal = !isLeaf;
      const isActiveInternal = isInternal && activeNode === u;

      const mark = completedByNode[u]; // true | Method | undefined
      const isCompletedForCurrent =
        mark === true || (!!activeMethod && mark === activeMethod);

      const showPlus = isActiveInternal;
      const showWhiteCap = isInternal && (isActiveInternal || !isCompletedForCurrent);

      const rOuter = isLeaf ? 3.75 : 3.25;

      out.push({
        id: u,
        x,
        y,
        rOuter,
        color: liveColor,
        showWhiteCap,
        showPlus,
        isLeaf,
      });
    }
    return out;
  }, [
    core.nodes,
    core.deg,
    daylight.angle,
    cx,
    cy,
    radiusOfDepth,
    completedByNode,
    activeMethod,
    activeNode,
    liveColor,
    weightedDepth,
  ]);

  // Total rendered edge length
  const totalRenderedLength = React.useMemo(() => {
    let sum = 0;
    for (const [v, p] of core.parent.entries()) {
      if (p === null) continue;
      const av = daylight.angle.get(v) ?? 0;
      const ap = daylight.angle.get(p) ?? 0;
      const dv = (weightedDepth.get(v) ?? 0);
      const dp = (weightedDepth.get(p) ?? 0);
      const pv = polar(cx, cy, radiusOfDepth(dv), av);
      const pp = polar(cx, cy, radiusOfDepth(dp), ap);
      sum += Math.hypot(pv.x - pp.x, pv.y - pp.y);
    }
    return sum;
  }, [core.parent, daylight.angle, weightedDepth, cx, cy, radiusOfDepth]);

  // Build the center-scaling transform string
  const contentTransform = React.useMemo(
    () => `translate(${cx},${cy}) scale(${zoom}) translate(${-cx},${-cy})`,
    [cx, cy, zoom]
  );

  const wsvg = Math.max(0, size.w);
  const hsvg = Math.max(0, size.h);

  return (
    <div ref={ref} className={`w-full h-full relative ${className}`}>
      <svg
        width={wsvg}
        height={hsvg}
        viewBox={`0 0 ${wsvg} ${hsvg}`}
        role="img"
        aria-label="Daylight-adjusted unrooted tree"
      >
        {/* Apply vertical pixel shift outside the zoom transform */}
        <g transform={`translate(0, ${yShift})`}>
          <g transform={contentTransform}>
            {/* Links */}
            <g stroke="#334155" strokeOpacity={0.25} strokeWidth={1}>
              {links.map((L) => (
                <line key={L.key} x1={L.x1} y1={L.y1} x2={L.x2} y2={L.y2} />
              ))}
            </g>

            {/* Nodes */}
            <g>
              {nodesViz.map((n) => (
                <circle
                  key={`base-${n.id}`}
                  cx={n.x}
                  cy={n.y}
                  r={n.rOuter}
                  fill={n.color}
                  stroke={n.color}
                  strokeWidth={1.2}
                />
              ))}

              {/* White cap for internal nodes that were NOT visited in any pass */}
              {nodesViz.map((n) =>
                n.showWhiteCap ? (
                  <circle
                    key={`cap-${n.id}`}
                    cx={n.x}
                    cy={n.y}
                    r={Math.max(1.6, n.rOuter * 0.8)}
                    fill="#ffffff"
                    stroke="#ffffff"
                    strokeWidth={0.6}
                  />
                ) : null
              )}

              {/* Plus symbol only for the active internal node while incomplete */}
              {nodesViz.map((n) =>
                n.showPlus ? (
                  <g key={`plus-${n.id}`}>
                    <line
                      x1={n.x - n.rOuter * 0.6}
                      y1={n.y}
                      x2={n.x + n.rOuter * 0.6}
                      y2={n.y}
                      stroke={n.color}
                      strokeWidth={1.6}
                      strokeLinecap="round"
                    />
                    <line
                      x1={n.x}
                      y1={n.y - n.rOuter * 0.6}
                      x2={n.x}
                      y2={n.y + n.rOuter * 0.6}
                      stroke={n.color}
                      strokeWidth={1.6}
                      strokeLinecap="round"
                    />
                  </g>
                ) : null
              )}
            </g>

            {/* Labels (both leaves and internal nodes) */}
            {showLabels ? (
              <>
                {/* Leaf labels: pushed outward along the ray */}
                <g fontSize={12} fill="#000">
                  {nodesViz
                    .filter((n) => n.isLeaf)
                    .map((n) => {
                      const theta = daylight.angle.get(n.id) ?? 0;
                      const dv = weightedDepth.get(n.id) ?? 0;
                      const R = radiusOfDepth(dv);
                      const dirX = Math.cos(theta - PI / 2);
                      const dirY = Math.sin(theta - PI / 2);
                      const offset = n.rOuter + 12;
                      const tx = n.x + dirX * offset;
                      const ty = n.y + dirY * offset;
                      const anchor =
                        dirX > 0.15 ? "start" : dirX < -0.15 ? "end" : "middle";
                      return (
                        <text
                          key={`leaf-label-${n.id}`}
                          x={tx}
                          y={ty}
                          textAnchor={anchor}
                          dominantBaseline="middle"
                        >
                          {n.id}
                        </text>
                      );
                    })}
                </g>

                {/* Internal node labels: centered just above each node */}
                <g fontSize={11} fill="#000">
                  {nodesViz
                    .filter((n) => !n.isLeaf)
                    .map((n) => (
                      <text
                        key={`internal-label-${n.id}`}
                        x={n.x}
                        y={n.y - (n.rOuter + 6)}
                        textAnchor="middle"
                        dominantBaseline="ideographic"
                      >
                        {n.id}
                      </text>
                    ))}
                </g>
              </>
            ) : null}
          </g>
        </g>
      </svg>

      {/* Side panel with diagnostics (not scaled) */}
      {/* <div className="absolute top-2 right-2 bg-white/85 backdrop-blur-sm border border-black/10 shadow-sm rounded-xl p-3 max-w-[280px] text-[12px] text-black leading-tight">
        <div className="font-semibold mb-1">Tree diagnostics</div>
        <div className="mb-1">
          <span className="opacity-70">Daylight passes:</span>{" "}
          <span>{Math.max(1, Math.floor(daylightPasses))}</span>
        </div>
        <div className="font-medium mb-1">Daylight imbalance (max−min gap)</div>
        <ul className="space-y-0.5 max-h-56 overflow-auto pr-1">
          {daylightImbalance.length === 0 ? (
            <li className="opacity-60">No internal nodes</li>
          ) : (
            daylightImbalance.map(({ id, imbalanceRad }) => (
              <li key={id} className="flex justify-between gap-2">
                <span className="opacity-80">{id}</span>
                <span>{(imbalanceRad * 180 / Math.PI).toFixed(1)}°</span>
              </li>
            ))
          )}
        </ul>
      </div> */}
    </div>
  );
}
