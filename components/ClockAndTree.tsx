// components/ClockAndTree.tsx
"use client";

import React from "react";
import StainedGlassClock from "@/components/StainedGlassClock";
import DaylightTree from "@/components/DaylightTree";
import { useProgressStore } from "@/stores/progress";

type Props = {
  clockHeight?: number;
  treeHeight?: number;
  className?: string;
  /** How many global passes of the daylight algorithm to run (default 1) */
  daylightPasses?: number;
  /** Vertical pixel shift applied to the tree (not scaled by zoom) */
  yShift?: number;
};

type EdgeInput =
  | [string, string, number]
  | { source: string; target: string; length?: number | null };

function toTuple(e: EdgeInput): [string, string] {
  if (Array.isArray(e)) {
    const [u, v] = e;
    return [u, v];
  }
  return [e.source, e.target];
}

export default function ClockAndTree({
  clockHeight = 520,
  treeHeight = 340,
  className = "",
  daylightPasses = 0,
  yShift = -60,
}: Props) {
  // Progress state
  const edges = useProgressStore((s) => s.edges) as EdgeInput[] | undefined;
  const currentMethod = useProgressStore((s) => s.currentMethod);
  const completedByNode = useProgressStore((s) => s.completedByNode);
  const totalReps = useProgressStore((s) => s.totalReps);
  const currentRep = useProgressStore((s) => s.currentRep);

  // Clock (drives the active node index)
  const longHandMethod = useProgressStore((s) => s.longHandMethod);
  const longHandNodeIndex = useProgressStore((s) => s.longHandNodeIndex);
  const longHandNodeTotal = useProgressStore((s) => s.longHandNodeTotal);

  // Build a stable list of INTERNAL nodes (degree > 1), sorted by id
  const internalNodes = React.useMemo(() => {
    if (!edges || edges.length === 0) return [] as string[];
    const deg = new Map<string, number>();
    for (const e of edges) {
      const [u, v] = toTuple(e);
      deg.set(u, (deg.get(u) || 0) + 1);
      deg.set(v, (deg.get(v) || 0) + 1);
    }
    return Array.from(deg.entries())
      .filter(([_, d]) => d > 1)
      .map(([id]) => id)
      .sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  }, [edges]);

  // Pick the single active internal node from the index
  const activeNode = React.useMemo(() => {
    if (!internalNodes.length) return null;
    if (typeof longHandNodeIndex !== "number") return null;
    const size = internalNodes.length;
    const idx = ((longHandNodeIndex % size) + size) % size; // safe modulo
    return internalNodes[idx] ?? null;
  }, [internalNodes, longHandNodeIndex]);

  // Show tree until *all* repetitions are finished
  const showTree = !!edges && edges.length > 0 && currentRep < totalReps;

  return (
    <div className={`flex flex-col gap-3 ${className}`}>
      {/* Clock */}
      <div className="w-full rounded-2xl" style={{ height: clockHeight }}>
        <StainedGlassClock
          method={longHandMethod}
          rep={currentRep}
          totalReps={totalReps}
          className="p-3"
        />
      </div>

      {/* Tree */}
      <div
        className="w-full rounded-2xl overflow-hidden relative"
        style={{ height: treeHeight }}
      >
        {showTree ? (
          <DaylightTree
            className="p-3"
            edges={edges || []}
            activeMethod={currentMethod}
            completedByNode={completedByNode}
            activeNode={activeNode}
            padding={14}
            showLabels={true}
            daylightPasses={daylightPasses}
            yShift={yShift}
          />
        ) : (
          <div className="absolute inset-0 grid place-items-center text-sm text-black/60" />
        )}
      </div>
    </div>
  );
}
