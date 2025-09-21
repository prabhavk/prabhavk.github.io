// components/StainedGlassClock.tsx
"use client";

import React from "react";
import type { Method } from "@/stores/progress";
import { useProgressStore } from "@/stores/progress";

type Props = {
  method: Method | null;       // optional override for long-hand method
  rep: number;                  // (unused for short hand; store-driven)
  totalReps: number;            // (fallback only; store-driven preferred)
  minuteIndex?: number;
  minuteTotal?: number;
  className?: string;
};

const TAU = Math.PI * 2;

export default function StainedGlassClock({
  method,
  rep,
  totalReps,
  className = "",
}: Props) {
  const ref = React.useRef<HTMLDivElement>(null);
  const [size, setSize] = React.useState({ w: 0, h: 0 });
  const markerId = React.useId();

  // --- Store state ---
  const colors = useProgressStore((s) => s.colors);

  // SHORT HAND: read progress directly from the store so it moves on every "rep"
  const storeRep        = useProgressStore((s) => s.currentRep);
  const storeTotalReps  = useProgressStore((s) => s.totalReps);
  const storeDefaultReps= useProgressStore((s) => s.emInput.reps);

  // LONG HAND (sticky, moves only on layer 2 completion)
  const longHandMethod     = useProgressStore((s) => s.longHandMethod);
  const longHandNodeIndex  = useProgressStore((s) => s.longHandNodeIndex);
  const longHandNodeTotal  = useProgressStore((s) => s.longHandNodeTotal);

  const sticky = React.useMemo(
    () => ({
      method: longHandMethod as Method | null,
      nodeIndex: longHandNodeIndex as number | null,
      nodeTotal: longHandNodeTotal as number | null,
    }),
    [longHandMethod, longHandNodeIndex, longHandNodeTotal]
  );

  const storeMethod = useProgressStore((s) => s.currentMethod);
  const edges = useProgressStore((s) => s.edges);

  React.useEffect(() => {
    if (!ref.current) return;
    const ro = new ResizeObserver((entries) => {
      for (const e of entries) setSize({ w: e.contentRect.width, h: e.contentRect.height });
    });
    ro.observe(ref.current);
    return () => ro.disconnect();
  }, []);

  const cx = size.w / 2;
  const cy = size.h / 2;
  const R = Math.max(0, Math.min(size.w, size.h) / 2 - 10);

  /* ---------- Internal node label helpers ---------- */

  const isInternalLike = (name: string) =>
    /^h[_-]?\d+$/i.test(name) || /^h_/i.test(name);

  // Force internal-like labels to lowercase 'h#'
  const prettyLabel = (name: string): string => {
    const m = String(name).match(/^[Hh][_-]?(\d+)$/);
    if (m) return `h${m[1]}`;          // <-- always lowercase h
    if (!name) return name;
    return name.charAt(0).toUpperCase() + name.slice(1);
  };

  const rootLabels = React.useMemo(() => {
    const set = new Set<string>();
    for (const e of edges || []) {
      const u = e?.[0]; const v = e?.[1];
      if (typeof u === "string" && isInternalLike(u)) set.add(u);
      if (typeof v === "string" && isInternalLike(v)) set.add(v);
    }
    if (set.size === 0) return ["h1", "h2", "h3", "h4", "h5", "h6"];
    const arr = Array.from(set).sort((a, b) => {
      const ma = a.match(/^[Hh][_-]?(\d+)$/);
      const mb = b.match(/^[Hh][_-]?(\d+)$/);
      if (ma && mb) return Number(ma[1]) - Number(mb[1]);
      return a.localeCompare(b);
    });
    return arr.map(prettyLabel);
  }, [edges]);

  /* ---------- Geometry ---------- */

  const startCenter = -TAU / 4;      // 12 o’clock
  const methodWedge = TAU / 3;
  const METHOD_ORDER: Method[] = ["dirichlet", "parsimony", "hss"];

  // Prefer sticky totals only when sticky exists; else derive from labels
  const hasSticky = sticky.method != null && sticky.nodeIndex != null;
  const N = Math.max(1, hasSticky && sticky.nodeTotal ? sticky.nodeTotal : rootLabels.length);
  const sliceStep = methodWedge / N;

  const sectorPath = React.useCallback(
    (a0: number, a1: number, r0: number, r1: number) => {
      const p = (ang: number, r: number) => [Math.cos(ang) * r, Math.sin(ang) * r] as const;
      const [x0o, y0o] = p(a0, r1);
      const [x1o, y1o] = p(a1, r1);
      const large = a1 - a0 > Math.PI ? 1 : 0;
      if (r0 <= 0.5) {
        return `M ${x0o} ${y0o} A ${r1} ${r1} 0 ${large} 1 ${x1o} ${y1o} L 0 0 Z`;
      }
      const [x1i, y1i] = p(a1, r0);
      const [x0i, y0i] = p(a0, r0);
      return [
        `M ${x0o} ${y0o}`,
        `A ${r1} ${r1} 0 ${large} 1 ${x1o} ${y1o}`,
        `L ${x1i} ${y1i}`,
        `A ${r0} ${r0} 0 ${large} 0 ${x0i} ${y0i}`,
        "Z",
      ].join(" ");
    },
    []
  );

  /* ---------- Hand angles ---------- */

  // SHORT HAND — drive from store so it only moves at end-of-rep:
  const denomReps =
    (storeTotalReps && storeTotalReps > 0 ? storeTotalReps :
     (totalReps && totalReps > 0 ? totalReps : (storeDefaultReps || 30)));

  const clampedRep = Math.max(0, Math.min(denomReps, Number.isFinite(storeRep) ? storeRep : 0));
  const shortAngle =
    clampedRep > 0 && denomReps > 0
      ? startCenter + (clampedRep / denomReps) * TAU
      : startCenter;

  // LONG HAND — end of slice for layer-3 completion (sticky)
  const methodUsed: Method =
    (hasSticky ? sticky.method : (method ?? storeMethod ?? "dirichlet")) as Method;
  const mi = Math.max(0, METHOD_ORDER.indexOf(methodUsed));
  const wedgeStart = startCenter + mi * methodWedge;

  const BOUNDARY_EPS = 1e-6;
  let longAngle = startCenter;
  if (hasSticky) {
    const idx0 = Math.min(N, Math.max(1, sticky.nodeIndex!)) - 1; // 0-based
    longAngle = wedgeStart + (idx0 + 1) * sliceStep - BOUNDARY_EPS; // end of slice
  }

  /* ---------- Hand geometry & styling ---------- */

  const baseShortLen = R * 0.50;
  const baseLongLen  = R * 0.86;

  const HEAD_LEN = Math.max(7, Math.min(14, R * 0.08));
  const HEAD_W   = Math.max(6, Math.min(12, R * 0.06));
  const OVERSHOOT = Math.max(2, Math.min(6, R * 0.03));

  const shortTipR = baseShortLen + OVERSHOOT;
  const longTipR  = baseLongLen  + OVERSHOOT;

  const shortLineR = Math.max(0, shortTipR - HEAD_LEN);
  const longLineR  = Math.max(0, longTipR  - HEAD_LEN);

  const shortEnd = { x: Math.cos(shortAngle) * shortLineR, y: Math.sin(shortAngle) * shortLineR };
  const longEnd  = { x: Math.cos(longAngle)  * longLineR,  y: Math.sin(longAngle)  * longLineR  };

  const SLICE_OPACITY = 0.85;
  const HAIRLINE = "rgba(0,0,0,0.06)";
  const TICK_COLOR = "rgba(0,0,0,0.35)";

  const labelOffset = Math.max(12, Math.min(28, R * 0.14));
  const viewPad = labelOffset + 10;

  const boundary = (angle: number, label: string) => {
    const r1 = R * 0.99;
    const rx = Math.cos(angle) * r1;
    const ry = Math.sin(angle) * r1;
    const lr = R + labelOffset * 0.6;
    const lx = Math.cos(angle) * lr;
    const ly = Math.sin(angle) * lr;
    return (
      <g key={`boundary-${label}`}>
        <line x1={0} y1={0} x2={rx} y2={ry} stroke={TICK_COLOR} strokeWidth={1} />
        <text
          x={lx}
          y={ly}
          textAnchor="middle"
          dominantBaseline="middle"
          fontSize={14}
          fill="#111827"
          style={{ fontWeight: 700 }}
        >
          {label}
        </text>
      </g>
    );
  };

  return (
    <div ref={ref} className={`w-full h-full ${className}`}>
      <svg
        width={size.w}
        height={size.h}
        viewBox={`${-viewPad} ${-viewPad} ${size.w + 2 * viewPad} ${size.h + 2 * viewPad}`}
      >
        <defs>
          <marker
            id={`arrow-head-black-${markerId}`}
            markerUnits="userSpaceOnUse"
            viewBox={`0 0 ${HEAD_LEN} ${HEAD_W}`}
            refX={HEAD_LEN}
            refY={HEAD_W / 2}
            markerWidth={HEAD_LEN}
            markerHeight={HEAD_W}
            orient="auto"
          >
            <path d={`M 0 0 L ${HEAD_LEN} ${HEAD_W / 2} L 0 ${HEAD_W} z`} fill="#000000" />
          </marker>
        </defs>

        <g transform={`translate(${cx},${cy})`}>
          {(["dirichlet", "parsimony", "hss"] as const).map((m, idx) => {
            const wedgeStartLocal = startCenter + idx * methodWedge;
            const fill = colors[m].middle;
            return (
              <g key={`wedge-${m}`}>
                {Array.from({ length: N }, (_, j) => {
                  const a0 = wedgeStartLocal + j * sliceStep;
                  const a1 = a0 + sliceStep;
                  const path = sectorPath(a0, a1, 0, R);

                  const ac = a0 + sliceStep / 2;
                  const tickR0 = R * 0.92;
                  const tickR1 = R * 0.99;
                  const x0 = Math.cos(ac) * tickR0, y0 = Math.sin(ac) * tickR0;
                  const x1 = Math.cos(ac) * tickR1, y1 = Math.sin(ac) * tickR1;

                  const lr = R + labelOffset * 0.45;
                  const lx = Math.cos(ac) * lr;
                  const ly = Math.sin(ac) * lr;

                  const label = rootLabels[j] ?? `h${j + 1}`; // <-- fallback lowercase

                  return (
                    <g key={`${m}-slice-${j}`}>
                      <path d={path} fill={fill} opacity={SLICE_OPACITY} stroke={HAIRLINE} strokeWidth={1} />
                      <line x1={x0} y1={y0} x2={x1} y2={y1} stroke={TICK_COLOR} strokeWidth={1} />
                      <text
                        x={lx}
                        y={ly}
                        textAnchor="middle"
                        dominantBaseline="middle"
                        fontSize={12}
                        fill={"#111827"}
                        opacity={0.9}
                        style={{ fontWeight: 600 }}
                      >
                        {label}
                      </text>
                    </g>
                  );
                })}
              </g>
            );
          })}

          {/* Method boundaries (D | P | H@12) */}
          {boundary(startCenter + 1 * methodWedge, "D")}
          {boundary(startCenter + 2 * methodWedge, "P")}
          {boundary(startCenter + 3 * methodWedge, "H")}

          {/* SHORT HAND (store-driven) */}
          <line
            x1={0}
            y1={0}
            x2={shortEnd.x}
            y2={shortEnd.y}
            stroke="#000000"
            strokeWidth={3}
            strokeLinecap="round"
            markerEnd={`url(#arrow-head-black-${markerId})`}
            opacity={0.95}
          />

          {/* LONG HAND (sticky: end of slice on layer-3 completion) */}
          <line
            x1={0}
            y1={0}
            x2={longEnd.x}
            y2={longEnd.y}
            stroke="#000000"
            strokeWidth={3}
            strokeLinecap="round"
            markerEnd={`url(#arrow-head-black-${markerId})`}
            opacity={0.95}
          />

          <circle cx={0} cy={0} r={Math.max(2, R * 0.02)} fill="#000" opacity={0.9} />
        </g>
      </svg>
    </div>
  );
}
