// components/tree-viewer.tsx
"use client";

import React, { useMemo } from "react";

/** ---- Minimal Newick parser (no deps) ---- */
type NwNode = {
  name?: string;
  length?: number;
  children?: NwNode[];
};

function parseNewick(input: string): NwNode {
  const s = (input || "").trim();
  let i = 0;

  const ws = () => { while (i < s.length && /\s/.test(s[i]!)) i++; };

  const readName = (): string => {
    ws();
    const start = i;
    while (i < s.length && !/[,\(\):;\s]/.test(s[i]!)) i++;
    return s.slice(start, i);
  };

  const readNumber = (): number | undefined => {
    ws();
    const start = i;
    while (i < s.length && /[0-9eE\+\-\.]/.test(s[i]!)) i++;
    if (i === start) return undefined;
    const v = Number(s.slice(start, i));
    return Number.isFinite(v) ? v : undefined;
  };

  const parseSubtree = (): NwNode => {
    ws();
    const node: NwNode = {};
    if (s[i] === "(") {
      i++; // '('
      node.children = [];
      for (;;) {
        node.children.push(parseSubtree());
        ws();
        if (s[i] === ",") { i++; continue; }
        if (s[i] === ")") { i++; break; }
        if (i >= s.length) throw new Error("Unclosed '(' in Newick");
      }
      ws();
      const nm = readName();
      if (nm) node.name = nm;
    } else {
      const nm = readName();
      if (nm) node.name = nm;
    }
    ws();
    if (s[i] === ":") {
      i++;
      const len = readNumber();
      if (len !== undefined) node.length = len;
    }
    return node;
  };

  const root = parseSubtree();
  ws();
  if (s[i] === ";") i++;
  return root;
}

/** ---- Simple unrooted radial layout ----
 * Strategy:
 *  - Place leaves equally spaced on a circle.
 *  - Place each internal node at the centroid of its descendant leaves' positions.
 *  - Draw straight segments parent↔child.
 */
type Point = { x: number; y: number };
type Positioned = NwNode & { id: number; parent?: Positioned | null; pos?: Point };

function collectLeaves(n: NwNode, out: NwNode[] = []): NwNode[] {
  if (!n.children || n.children.length === 0) out.push(n);
  else for (const c of n.children) collectLeaves(c, out);
  return out;
}

function annotateParent(n: Positioned, parent: Positioned | null = null, nextId = { v: 0 }) {
  n.id = nextId.v++;
  n.parent = parent;
  if (n.children) {
    for (const c of n.children as Positioned[]) annotateParent(c, n, nextId);
  }
}

function leavesOf(n: NwNode): NwNode[] {
  return collectLeaves(n, []);
}

function positionUnrooted(root: NwNode, width: number, height: number, pad = 20) {
  const W = Math.max(200, width);
  const H = Math.max(200, height);
  const cx = W / 2, cy = H / 2;
  const radius = Math.max(10, Math.min(W, H) / 2 - pad);

  // Clone shallowly so we can annotate
  const clone = (n: NwNode): Positioned => ({
    name: n.name,
    length: n.length,
    children: n.children ? n.children.map(clone) : undefined,
  }) as Positioned;
  const R = clone(root);
  annotateParent(R, null);

  // 1) place leaves evenly around circle
  const leaves = leavesOf(R);
  const nLeaves = Math.max(1, leaves.length);
  const angleForLeaf = new Map<NwNode, number>();
  for (let k = 0; k < nLeaves; k++) {
    const a = (2 * Math.PI * k) / nLeaves;
    angleForLeaf.set(leaves[k]!, a);
  }

  // 2) postorder walk to compute descendant leaf angles
  const leafAngleList = new Map<Positioned, number[]>();
  const post = (n: Positioned): number[] => {
    if (!n.children || n.children.length === 0) {
      const a = angleForLeaf.get(n as NwNode) ?? 0;
      leafAngleList.set(n, [a]);
      return [a];
    }
    const all: number[] = [];
    for (const c of n.children as Positioned[]) {
      all.push(...post(c));
    }
    leafAngleList.set(n, all);
    return all;
  };
  post(R);

  // 3) assign positions
  const posOf = new Map<Positioned, Point>();
  const setPos = (n: Positioned) => {
    const list = leafAngleList.get(n) || [];
    if (!n.children || n.children.length === 0) {
      const a = list[0] ?? 0;
      posOf.set(n, { x: cx + radius * Math.cos(a), y: cy + radius * Math.sin(a) });
      return;
    }
    // centroid of descendant leaves on circle
    let sx = 0, sy = 0;
    for (const a of list) { sx += Math.cos(a); sy += Math.sin(a); }
    const denom = Math.max(1e-9, list.length);
    const px = cx + radius * (sx / denom);
    const py = cy + radius * (sy / denom);
    posOf.set(n, { x: px, y: py });
    for (const c of (n.children || []) as Positioned[]) setPos(c);
  };
  setPos(R);

  // 4) edges parent<->child, and label positions for leaves
  const edges: Array<{ x1: number; y1: number; x2: number; y2: number }> = [];
  const labels: Array<{ x: number; y: number; text: string }> = [];
  const walk = (n: Positioned) => {
    const p = posOf.get(n)!;
    if (!n.children || n.children.length === 0) {
      labels.push({ x: p.x, y: p.y, text: n.name || "" });
    } else {
      for (const c of (n.children || []) as Positioned[]) {
        const q = posOf.get(c)!;
        edges.push({ x1: p.x, y1: p.y, x2: q.x, y2: q.y });
        walk(c);
      }
    }
  };
  walk(R);

  return { width: W, height: H, edges, labels };
}

/** ---- Component ---- */
type Props = {
  newick?: string | null;
  height?: number | string; // e.g. 520 or "70vh"
  className?: string;
  // optional fixed width (defaults to 900 if you pass a number height)
  width?: number;
};

export default function TreeViewer({ newick, height = 520, className, width }: Props) {
  // Normalize Newick once
  const normalized = useMemo(() => {
    if (typeof newick !== "string") return "";
    const t = newick.trim();
    return t && !t.endsWith(";") ? `${t};` : t;
  }, [newick]);

  // Parse (no throwing outside)
  const { root, error } = useMemo(() => {
    if (!normalized) return { root: null as NwNode | null, error: "No Newick string provided." };
    try {
      return { root: parseNewick(normalized), error: null as string | null };
    } catch (e) {
      const msg = e instanceof Error ? e.message : String(e);
      return { root: null, error: `Failed to parse Newick: ${msg}` };
    }
  }, [normalized]);

  // Compute layout
  const layout = useMemo(() => {
    if (!root) return null;
    const H = typeof height === "number" ? height : 520;
    const W = typeof width === "number" ? width : 900;
    return positionUnrooted(root, W, H);
  }, [root, height, width]);

  // Render
  const style: React.CSSProperties = {
    width: "100%",
    height: typeof height === "number" ? `${height}px` : height,
    border: "1px solid #eee",
    background: "#fff",
    borderRadius: 4,
  };

  if (error) {
    return <div style={style} className={className}><div style={{ padding: 12, color: "crimson" }}>{error}</div></div>;
  }
  if (!layout) {
    return <div style={style} className={className}><div style={{ padding: 12 }}>No tree to display.</div></div>;
  }

  const { width: W, height: H, edges, labels } = layout;

  return (
    <div style={style} className={className}>
      {/* viewBox scales to container width; edges are precomputed in absolute coords */}
      <svg viewBox={`0 0 ${W} ${H}`} width="100%" height="100%">
        <g stroke="#111" strokeWidth={1}>
          {edges.map((e, i) => (
            <line key={i} x1={e.x1} y1={e.y1} x2={e.x2} y2={e.y2} />
          ))}
        </g>
        <g fill="#111" fontSize={12}>
          {labels.map((l, i) => (
            <text key={i} x={l.x + 6} y={l.y + 4} style={{ userSelect: "none" }}>
              {l.text}
            </text>
          ))}
        </g>
      </svg>
    </div>
  );
}
