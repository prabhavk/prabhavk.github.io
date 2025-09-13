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
  showLabels?: boolean;
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
};
function rootAtNode(adj: Adj, rootId: string): Rooted {
  if (!rootId) return { root: "", nodes: [], parent: new Map(), children: new Map(), deg: new Map() };

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
  return { root: rootId, nodes: Array.from(children.keys()), parent, children, deg };
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
  // collect leaves in DFS order
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

/* ---------- Euler tour helpers (for subtree membership) ---------- */
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

/* ================================================================
   =======  PUBLIC HELPERS YOU ASKED FOR (exported)  ==============
   ================================================================ */

/** Wedge of a subtree adjacent to node p, described by left/right rays (angles at p). */
export type SubtreeWedge = {
  /** "child" for a child subtree, "parent" for the complement of subtree(p) */
  kind: "child" | "parent";
  /** id of the child root if kind === "child" */
  id?: string;
  /** leftmost and rightmost rays from p to leaves of this component (angles in [0, 2π)) */
  left: number;
  right: number;
  /** derived convenience values */
  center: number;
  width: number;
};

/**
 * Compute wedges (left/right rays) for all components attached to node `p`.
 * - Uses **equal-angle** node positions (passed via `angle` and `depth` + `radiusOfDepth`).
 * - Children components are each child subtree of `p`.
 * - If `p` has a parent, includes one "parent" wedge for all leaves **outside** subtree(p).
 *
 * NOTE: This function requires Euler tin/tout so it can test subtree membership.
 */
export function computeSubtreeWedgesAtNode(
  p: string,
  params: {
    parent: Map<string, string | null>;
    children: Map<string, string[]>;
    leaves: string[];                         // leaves of the whole tree
    angle: Map<string, number>;               // equal-angle angles for nodes
    depth: Map<string, number>;               // topological depths
    radiusOfDepth: (d: number) => number;     // converts depth -> radial distance
    cx: number;
    cy: number;
    tin: Map<string, number>;
    tout: Map<string, number>;
  }
): SubtreeWedge[] {
  const { parent, children, leaves, angle, depth, radiusOfDepth, cx, cy, tin, tout } = params;

  // position helper under equal-angle geometry
  const nodePos = (u: string) => {
    const a = angle.get(u) ?? 0;
    const d = depth.get(u) ?? 0;
    return polar(cx, cy, radiusOfDepth(d), a);
  };

  // direction (angle) of ray from p to q
  const rayAngleFromP = (q: string) => {
    const P = nodePos(p);
    const Q = nodePos(q);
    return as01Angle(Math.atan2(Q.y - P.y, Q.x - P.x));
  };

  // Collect neighbor components: each child, and the "parent side" if any
  const nbrs: Array<{ kind: "child" | "parent"; id?: string; dir: number }> = [];
  for (const c of (children.get(p) || [])) {
    nbrs.push({ kind: "child", id: c, dir: rayAngleFromP(c) });
  }
  const par = parent.get(p);
  if (par !== null && par !== undefined) {
    nbrs.push({ kind: "parent", dir: rayAngleFromP(par!) });
  }

  // Sort neighbors by the direction of their connecting edge (cyclic order around p)
  nbrs.sort((a, b) => a.dir - b.dir);

  // For each neighbor component, gather the angles of rays from p to its leaves.
  const wedges: SubtreeWedge[] = [];
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
      // leaves inside child subtree
      const ls = leaves.filter((L) => isInSubtree(tin, tout, L, nb.id!));
      leafAngles = raysToLeaves(ls);
    } else {
      // "parent" = all leaves outside subtree(p)
      const ls = leaves.filter((L) => !isInSubtree(tin, tout, L, p));
      leafAngles = raysToLeaves(ls);
    }

    // If none (degenerate), fall back to the neighbor edge direction
    if (!leafAngles.length) leafAngles = [nb.dir];

    // Shift each angle to lie near the neighbor direction (prevents 0/2π wrap from splitting the wedge)
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

/** Minimal circular distance between two angles (in [0, 2π)), result in [0, π]. */
function angleDistance(a: number, b: number): number {
  const d = Math.abs(as01Angle(a) - as01Angle(b));
  return Math.min(d, TAU - d);
}

/**
 * Daylight angle between two wedges attached to the same node.
 * Defined as the smaller of:
 *   1) distance(left_A, right_B)
 *   2) distance(right_A, left_B)
 */
export function daylightAngleBetween(A: SubtreeWedge, B: SubtreeWedge): number {
  const d1 = angleDistance(A.left, B.right);
  const d2 = angleDistance(A.right, B.left);
  return Math.min(d1, d2);
}

/* ================================================================
   =====================  COMPONENT (EA ONLY)  =====================
   ================================================================ */

export default function DaylightTree({
  className = "",
  edges,
  activeMethod = null,
  completedByNode,
  activeNode = null,
  padding = 12,
  showLabels = false,
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
    };
  }, [edges]);

  /* ---------- Equal-angle starting angles (no daylight refinement) ---------- */
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
      d === 0 ? 0 : ringInner + (ringOuter - ringInner) * (d / Math.max(1, equalAngle.maxDepth)),
    [ringInner, ringOuter, equalAngle.maxDepth]
  );

  // Live color always follows the active method (or neutral).
  const liveColor = React.useMemo(
    () => (activeMethod ? METHOD_COLOR[activeMethod] : NEUTRAL),
    [activeMethod]
  );

  // Links
  const links = React.useMemo(() => {
    const out: Array<{ x1: number; y1: number; x2: number; y2: number; key: string }> = [];
    for (const [v, p] of core.parent.entries()) {
      if (p === null) continue;
      const av = equalAngle.angle.get(v) ?? 0;
      const ap = equalAngle.angle.get(p) ?? 0;
      const dv = equalAngle.depth.get(v) ?? 0;
      const dp = equalAngle.depth.get(p) ?? 0;
      const pv = polar(cx, cy, radiusOfDepth(dv), av);
      const pp = polar(cx, cy, radiusOfDepth(dp), ap);
      out.push({ x1: pv.x, y1: pv.y, x2: pp.x, y2: pp.y, key: `${p}->${v}` });
    }
    return out;
  }, [core.parent, equalAngle.angle, equalAngle.depth, cx, cy, radiusOfDepth]);

  // Nodes: base colored circle for all; leaves solid; internal nodes white-cap until completed;
  // show "+" only on the active internal node while incomplete.
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
      const a = equalAngle.angle.get(u) ?? 0;
      const d = equalAngle.depth.get(u) ?? 0;
      const { x, y } = polar(cx, cy, radiusOfDepth(d), a);

      const isLeaf = (core.deg.get(u) || 0) <= 1;
      const isActiveInternal = !isLeaf && activeNode === u;

      const mark = completedByNode[u]; // true | Method | undefined
      const markMethod = typeof mark === "string" ? (mark as Method) : null;

      const isCompleteForCurrent =
        mark === true || (!!activeMethod && markMethod === activeMethod && !isActiveInternal);

      const rOuter = isLeaf ? 3.75 : 3.25;

      out.push({
        id: u,
        x,
        y,
        rOuter,
        color: liveColor,
        showWhiteCap: !isLeaf && !isCompleteForCurrent,
        showPlus: isActiveInternal && !isCompleteForCurrent,
        isLeaf,
      });
    }
    return out;
  }, [
    core.nodes,
    core.deg,
    equalAngle.angle,
    equalAngle.depth,
    cx,
    cy,
    radiusOfDepth,
    completedByNode,
    activeMethod,
    activeNode,
    liveColor,
  ]);

  const wsvg = Math.max(0, size.w);
  const hsvg = Math.max(0, size.h);

  return (
    <div ref={ref} className={`w-full h-full ${className}`}>
      <svg
        width={wsvg}
        height={hsvg}
        viewBox={`0 0 ${wsvg} ${hsvg}`}
        role="img"
        aria-label="Equal-angle unrooted tree"
      >
        {/* Links */}
        <g stroke="#334155" strokeOpacity={0.25} strokeWidth={1}>
          {links.map((L) => (
            <line key={L.key} x1={L.x1} y1={L.y1} x2={L.x2} y2={L.y2} />
          ))}
        </g>

        {/* Nodes: base solid for all */}
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

          {/* White cap only for internal, incomplete nodes */}
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
                  stroke={liveColor}
                  strokeWidth={1.6}
                  strokeLinecap="round"
                />
                <line
                  x1={n.x}
                  y1={n.y - n.rOuter * 0.6}
                  x2={n.x}
                  y2={n.y + n.rOuter * 0.6}
                  stroke={liveColor}
                  strokeWidth={1.6}
                  strokeLinecap="round"
                />
              </g>
            ) : null
          )}
        </g>

        {/* Labels (optional) */}
        {showLabels ? (
          <g fontSize={11} fill="#111827" fillOpacity={0.85}>
            {nodesViz.map((n) => (
              <text
                key={`label-${n.id}`}
                x={n.x}
                y={n.y - (n.rOuter + 4)}
                textAnchor="middle"
                dominantBaseline="ideographic"
              >
                {n.id}
              </text>
            ))}
          </g>
        ) : null}
      </svg>
    </div>
  );
}
