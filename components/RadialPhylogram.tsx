// components/RadialPhylogram.tsx
"use client";

import * as d3 from "d3";
import React, { useEffect, useMemo, useRef } from "react";

export type Edge = { source: string; target: string; length?: number | null };

type Props = {
  edges: Edge[];
  root: string;
  width?: number;
  height?: number;
  padding?: number;
  showEdgeLengths?: boolean;
  daylightGapDeg?: number;
};

type NodeDatum = { name: string; bl?: number; children?: NodeDatum[] };
type H = d3.HierarchyNode<NodeDatum> & {
  sumLen?: number;
  tipCount?: number;
  theta0?: number;
  theta1?: number;
};
// same node but with the x/y we assign for linkRadial
type HPoint = H & { x: number; y: number };
type HLink = { source: HPoint; target: HPoint };
type MidLabel = { angle: number; radius: number; text: string };

// ---- build a rooted tree object from an undirected edge list ----
function buildTreeFromEdges(edges: Edge[], root: string): NodeDatum {
  const adj = new Map<string, Array<{ n: string; bl: number }>>();
  const add = (a: string, b: string, bl: number) => {
    if (!adj.has(a)) adj.set(a, []);
    adj.get(a)!.push({ n: b, bl });
  };
  for (const e of edges) {
    const a = String(e.source);
    const b = String(e.target);
    const bl = Number(e.length ?? 0) || 0;
    add(a, b, bl);
    add(b, a, bl);
  }
  if (!adj.has(root)) throw new Error(`Root "${root}" not present in edge list`);

  const seen = new Set<string>();
  function build(name: string): NodeDatum {
    seen.add(name);
    const node: NodeDatum = { name, bl: 0, children: [] };
    for (const { n, bl } of adj.get(name) || []) {
      if (!seen.has(n)) {
        const child = build(n);
        child.bl = bl; // branch length to parent
        node.children!.push(child);
      }
    }
    if (!node.children!.length) delete node.children;
    return node;
  }
  return build(root);
}

// ---- equal-daylight radial layout (angles), phylogram radius (lengths) ----
function layoutEqualDaylight(hroot: H, R: number, daylightGapRad: number) {
  // 1) tip counts
  const tips = (d: H): number => {
    if (!d.children || d.children.length === 0) return (d.tipCount = 1);
    d.tipCount = d.children.reduce((s, c) => s + tips(c as H), 0);
    return d.tipCount;
  };
  tips(hroot);

  // 2) cumulative branch length -> radius scale
  let maxLen = 0;
  hroot.eachBefore((d: H) => {
    const bl = Number(d.data.bl || 0);
    d.sumLen = (d.parent ? (d.parent as H).sumLen! : 0) + bl;
    if (d.sumLen! > maxLen) maxLen = d.sumLen!;
  });
  const rScale = d3.scaleLinear().domain([0, maxLen || 1]).range([0, R]);

  // 3) assign angular wedges [theta0, theta1] per node (equal-daylight)
  const assign = (d: H, a0: number, a1: number) => {
    d.theta0 = a0;
    d.theta1 = a1;
    // angle for this node (midpoint)
    (d as HPoint).x = (a0 + a1) / 2;
    // radius for this node from length
    (d as HPoint).y = rScale(d.sumLen || 0);

    if (!d.children || d.children.length === 0) return;

    const kids = d.children as H[];
    const totalTips = d.tipCount || kids.reduce((s, c) => s + (c.tipCount || 1), 0);

    const nGaps = Math.max(0, kids.length - 1);
    const daylightTotal = daylightGapRad * nGaps;
    const span = Math.max(0, (a1 - a0) - daylightTotal);

    // stable ordering; customize if needed
    kids.sort((a, b) => (a.data.name < b.data.name ? -1 : 1));

    let cursor = a0;
    for (let i = 0; i < kids.length; i++) {
      const c = kids[i];
      const frac = (c.tipCount || 1) / (totalTips || 1);
      const w = span * frac;
      const c0 = cursor;
      const c1 = cursor + w;
      assign(c, c0, c1);
      cursor = c1 + (i < kids.length - 1 ? daylightGapRad : 0);
    }
  };
  assign(hroot, 0, 2 * Math.PI);

  return { rScale };
}

export default function RadialPhylogram({
  edges,
  root,
  width = 700,
  height = 700,
  padding = 80,
  showEdgeLengths = false,
  daylightGapDeg = 1.5,
}: Props) {
  const ref = useRef<SVGSVGElement | null>(null);

  // Build hierarchy once per input
  const hroot = useMemo(() => {
    const treeObj = buildTreeFromEdges(edges, root);
    return d3.hierarchy<NodeDatum>(treeObj, (d) => d.children) as H;
  }, [edges, root]);

  useEffect(() => {
    if (!ref.current) return;
    const svg = d3.select(ref.current);
    svg.selectAll("*").remove();

    const W = width;
    const Ht = height;
    const R = Math.min(W, Ht) / 2 - padding;
    const daylightGapRad = (Math.PI / 180) * Math.max(0, daylightGapDeg);

    // equal-daylight for angle; length-scaled radius
    const { rScale } = layoutEqualDaylight(hroot, R, daylightGapRad);

    const g = svg
      .attr("viewBox", [-W / 2, -Ht / 2, W, Ht].join(" "))
      .attr("width", W)
      .attr("height", Ht)
      .append("g");

    // Build typed links with {x,y} present
    const links: HLink[] = hroot.links().map((l) => ({
      source: l.source as HPoint,
      target: l.target as HPoint,
    }));

    // Links
    const linkGen = d3
      .linkRadial<HLink, HPoint>()
      .angle((d) => d.x)
      .radius((d) => d.y);

    g.append("g")
      .selectAll<SVGPathElement, HLink>("path")
      .data(links)
      .join("path")
      .attr("d", (d) => linkGen(d) ?? "")
      .attr("fill", "none")
      .attr("stroke", "#111")
      .attr("stroke-width", 1.25);

    // Nodes
    const nodesSel = g
      .append("g")
      .selectAll<SVGGElement, HPoint>("g.node")
      .data(hroot.descendants().map((d) => d as HPoint))
      .join("g")
      .attr("class", "node")
      .attr(
        "transform",
        (d) => `rotate(${(d.x * 180) / Math.PI - 90}) translate(${d.y},0)`
      );

    nodesSel.append("circle").attr("r", 2.2).attr("fill", "#111");

    // Labels
    nodesSel
      .append("text")
      .attr("dy", "0.32em")
      .attr("x", (d) => (d.x < Math.PI ? 6 : -6))
      .attr("text-anchor", (d) => (d.x < Math.PI ? "start" : "end"))
      .attr("transform", (d) => (d.x >= Math.PI ? "rotate(180)" : null))
      .style(
        "font",
        "12px ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial"
      )
      .text((d) => d.data.name ?? "");

    // Edge-length labels (midpoint)
    if (showEdgeLengths) {
      const mids: MidLabel[] = links.map((l) => {
        const a = l.source;
        const b = l.target;
        return {
          angle: (a.x + b.x) / 2,
          radius: (a.y + b.y) / 2,
          text: (b.data?.bl ?? 0).toFixed(4),
        };
      });

      g.append("g")
        .selectAll<SVGTextElement, MidLabel>("text.edge-label")
        .data(mids)
        .join("text")
        .attr("class", "edge-label")
        .attr(
          "transform",
          (d) =>
            `rotate(${(d.angle * 180) / Math.PI - 90}) translate(${d.radius},0)`
        )
        .attr("dy", "-0.3em")
        .attr("text-anchor", "middle")
        .style("font", "10px ui-sans-serif, system-ui")
        .style("fill", "#555")
        .text((d) => d.text);
    }

    // Optional radial axis (branch length scale)
    const axis = d3.axisRight(rScale.copy().range([0, R])).ticks(5);
    svg
      .append("g")
      .attr("transform", `translate(${W / 2 - 20}, ${-Ht / 2 + 20})`)
      .call(axis)
      .call((gg) => gg.selectAll("text").style("font", "10px ui-sans-serif"));
  }, [hroot, width, height, padding, showEdgeLengths, daylightGapDeg]);

  return <svg ref={ref} />;
}
