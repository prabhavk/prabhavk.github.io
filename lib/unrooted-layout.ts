// lib/unrooted-layout.ts

/* ---------- Newick parser (minimal, no deps) ---------- */

export type NwNode = { name?: string; length?: number; children?: NwNode[] };

export function parseNewick(input: string): NwNode {
  const s = (input || "").trim();
  if (!s) throw new Error("Empty Newick");
  let i = 0;

  const isWS = (c: string) => /\s/.test(c);
  const skipWs = () => { while (i < s.length && isWS(s[i]!)) i++; };

  const readName = (): string => {
    skipWs();
    const start = i;
    while (i < s.length && !/[,\(\):;\s]/.test(s[i]!)) i++;
    return s.slice(start, i);
  };

  const readNumber = (): number | undefined => {
    skipWs();
    const start = i;
    while (i < s.length && /[0-9eE+\-.]/.test(s[i]!)) i++;
    if (i === start) return undefined;
    const v = Number(s.slice(start, i));
    return Number.isFinite(v) ? v : undefined;
  };

  const subtree = (): NwNode => {
    skipWs();
    const node: NwNode = {};
    if (s[i] === "(") {
      i++; // '('
      node.children = [];
      for (;;) {
        node.children.push(subtree());
        skipWs();
        if (s[i] === ",") { i++; continue; }
        if (s[i] === ")") { i++; break; }
        if (i >= s.length) throw new Error("Unclosed '('");
      }
      // optional internal node name
      const nm = readName();
      if (nm) node.name = nm;
    } else {
      const nm = readName();
      if (nm) node.name = nm;
    }
    skipWs();
    if (s[i] === ":") {
      i++;
      const len = readNumber();
      if (len !== undefined) node.length = len;
    }
    return node;
  };

  const root = subtree();
  skipWs();
  if (s[i] === ";") i++;
  return root;
}

/* ---------- Graph + unrooted equal-angle layout ---------- */

type Edge = { u: number; v: number; w: number }; // undirected, weight=w
type NodeRec = { id: number; name?: string; isLeaf: boolean };

type Graph = {
  nodes: NodeRec[];
  adj: Map<number, Array<{ v: number; w: number }>>;
};

function buildGraph(tree: NwNode): Graph {
  const nodes: NodeRec[] = [];
  const adj = new Map<number, Array<{ v: number; w: number }>>();

  const addNode = (name?: string): number => {
    const id = nodes.length;
    nodes.push({ id, name, isLeaf: true });
    adj.set(id, []);
    return id;
  };
  const addEdge = (a: number, b: number, w: number) => {
    adj.get(a)!.push({ v: b, w });
    adj.get(b)!.push({ v: a, w });
  };

  function build(n: NwNode, parentId: number | null): number {
    const id = addNode(n.name);
    const lenToParent = n.length ?? 1;

    if (parentId !== null) addEdge(parentId, id, lenToParent);

    if (n.children && n.children.length) {
      nodes[id].isLeaf = false;
      for (const c of n.children) build(c, id);
    }
    return id;
  }

  build(tree, null);
  return { nodes, adj };
}

function dijkstra(adj: Graph["adj"], start: number) {
  const N = adj.size;
  const dist = new Array<number>(N).fill(Infinity);
  const prev = new Array<number | null>(N).fill(null);
  const used = new Array<boolean>(N).fill(false);
  dist[start] = 0;

  for (let iter = 0; iter < N; iter++) {
    let u = -1, best = Infinity;
    for (let i = 0; i < N; i++) if (!used[i] && dist[i] < best) { best = dist[i]; u = i; }
    if (u === -1) break;
    used[u] = true;
    const nbrs = adj.get(u)!;
    for (const { v, w } of nbrs) {
      const nd = dist[u] + (Number.isFinite(w) ? w : 1);
      if (nd < dist[v]) { dist[v] = nd; prev[v] = u; }
    }
  }
  return { dist, prev };
}

function reconstructPath(prev: Array<number | null>, end: number): number[] {
  const path: number[] = [];
  for (let v: number | null = end; v !== null; v = prev[v]!) path.push(v);
  path.reverse();
  return path;
}

function findDiameter(adj: Graph["adj"]) {
  // 1) farthest from 0
  const a0 = dijkstra(adj, 0);
  let s = 0; for (let i = 0; i < a0.dist.length; i++) if (a0.dist[i] > a0.dist[s]) s = i;
  // 2) farthest from s
  const a1 = dijkstra(adj, s);
  let t = s; for (let i = 0; i < a1.dist.length; i++) if (a1.dist[i] > a1.dist[t]) t = i;
  const path = reconstructPath(a1.prev, t);
  const length = a1.dist[t];
  return { s, t, path, length };
}

function splitEdgeForMidpoint(
  g: Graph,
  path: number[],
  pathDist: number,
): number {
  // place root at midpoint of diameter
  const mid = pathDist / 2;
  // walk along path distances to find edge crossing midpoint
  // We'll need distances along path with edge weights
  const edgeW = (a: number, b: number) =>
    g.adj.get(a)!.find((x) => x.v === b)!.w;

  let accum = 0;
  for (let i = 0; i < path.length - 1; i++) {
    const a = path[i], b = path[i + 1];
    const w = edgeW(a, b);
    if (accum + w >= mid) {
      // Midpoint falls on edge (a,b): split by inserting synthetic root R
      const R = g.nodes.length;
      g.nodes.push({ id: R, name: undefined, isLeaf: false });
      g.adj.set(R, []);

      const wA = mid - accum;
      const wB = w - wA;

      // remove a<->b
      const rm = (p: number, q: number) => {
        const list = g.adj.get(p)!;
        const idx = list.findIndex((e) => e.v === q);
        if (idx >= 0) list.splice(idx, 1);
      };
      rm(a, b); rm(b, a);

      // connect a<->R (wA) and b<->R (wB)
      g.adj.get(a)!.push({ v: R, w: wA });
      g.adj.get(R)!.push({ v: a, w: wA });
      g.adj.get(b)!.push({ v: R, w: wB });
      g.adj.get(R)!.push({ v: b, w: wB });

      return R;
    }
    accum += w;
  }
  // If we didn’t split (numerical edge cases), choose center node
  return path[(path.length / 2) | 0];
}

type Rooted = {
  parent: number | null;
  children: number[];
  distFromRoot: number; // cumulative length
  leafCount: number;
  angle: number; // filled later
};
type RootedMap = Map<number, Rooted>;

function buildRootedFrom(g: Graph, root: number): RootedMap {
  const rooted: RootedMap = new Map();
  const seen = new Set<number>();

  function dfs(u: number, p: number | null, cum: number) {
    seen.add(u);
    const kids: number[] = [];
    for (const { v, w } of g.adj.get(u)!) {
      if (v === p) continue;
      kids.push(v);
      dfs(v, u, cum + w);
    }
    rooted.set(u, { parent: p, children: kids, distFromRoot: cum, leafCount: 0, angle: 0 });
  }
  dfs(root, null, 0);

  // leaf counts
  function count(u: number): number {
    const r = rooted.get(u)!;
    if (r.children.length === 0) {
      r.leafCount = 1;
      return 1;
    }
    let sum = 0;
    for (const v of r.children) sum += count(v);
    r.leafCount = sum;
    return sum;
  }
  count(root);
  return rooted;
}

function assignAnglesEqualAngle(root: number, rooted: RootedMap) {
  const TWO_PI = Math.PI * 2;
  function walk(u: number, start: number, span: number) {
    const r = rooted.get(u)!;
    if (r.children.length === 0) {
      r.angle = start + span / 2;
      return;
    }
    let cursor = start;
    for (const v of r.children) {
      const child = rooted.get(v)!;
      const frac = child.leafCount / r.leafCount;
      const childSpan = span * frac;
      walk(v, cursor, childSpan);
      cursor += childSpan;
    }
    // internal node angle: average of children angles
    const angles = r.children.map((v) => rooted.get(v)!.angle);
    r.angle = angles.reduce((a, b) => a + b, 0) / angles.length;
  }
  walk(root, 0, TWO_PI);
}

export type LayoutPoint = {
  id: number;
  name?: string;
  isLeaf: boolean;
  x: number;
  y: number;
  r: number;        // radius from center (scaled)
  angle: number;    // angle in radians
};
export type LayoutEdge = { u: number; v: number };
export type UnrootedLayout = {
  nodes: LayoutPoint[];
  edges: LayoutEdge[];
  center: { x: number; y: number };
  bounds: { width: number; height: number };
};

export function layoutUnrooted(
  newick: string,
  opts: { width?: number; height?: number; padding?: number } = {}
): UnrootedLayout {
  const width = Math.max(10, opts.width ?? 800);
  const height = Math.max(10, opts.height ?? 600);
  const pad = Math.max(0, opts.padding ?? 24);

  // 1) Parse, build undirected graph
  const t = newick && !newick.trim().endsWith(";") ? newick.trim() + ";" : (newick || "");
  const tree = parseNewick(t);
  const g = buildGraph(tree);

  // 2) Find diameter and root at its midpoint
  const diam = findDiameter(g.adj);
  const root = splitEdgeForMidpoint(g, diam.path, diam.length);

  // 3) Build rooted representation from that root
  const rooted = buildRootedFrom(g, root);

  // 4) Equal-angle assignment (unrooted style)
  assignAnglesEqualAngle(root, rooted);

  // 5) Radii: cumulative branch length from root -> scale to fit
  const maxDist = Math.max(...Array.from(rooted.values()).map((r) => r.distFromRoot));
  const R = Math.max(1e-9, Math.min((width - 2 * pad), (height - 2 * pad)) / 2);
  const scale = maxDist > 0 ? (R / maxDist) : 1;

  // 6) Coordinates
  const cx = width / 2, cy = height / 2;
  const nodes: LayoutPoint[] = g.nodes.map((n) => {
    const r = rooted.get(n.id)!;
    const rr = r.distFromRoot * scale;
    const x = cx + rr * Math.cos(r.angle);
    const y = cy + rr * Math.sin(r.angle);
    return { id: n.id, name: n.name, isLeaf: n.isLeaf, x, y, r: rr, angle: r.angle };
  });

  const edges: LayoutEdge[] = [];
  for (const [u, nbrs] of g.adj.entries()) {
    for (const { v } of nbrs) if (u < v) edges.push({ u, v });
  }

  return { nodes, edges, center: { x: cx, y: cy }, bounds: { width, height } };
}
