// components/CakeChart.tsx
"use client";
import React, { useMemo } from "react";

export type MethodName = "parsimony" | "dirichlet" | "hss";

type OneLayer = {
  layer: 0 | 1 | 2;     // 0=coarse, 1=medium, 2=fine
  iter: number | null;  // stripes count
  ll_final: number | null; // height driver
};

type LayersByMethod = Record<MethodName, OneLayer[]>;

type PatternCounts = { coarse: number; medium: number; fine: number };

type Props = {
  data: LayersByMethod;           // 3 entries: parsimony, dirichlet, hss; each with up to 3 layers
  patternCounts: PatternCounts;   // cumulative %
  rep: number | null;             // for labeling only
  root: string | null;            // for labeling only
};

// --- Shared color scheme (method -> color) ---
const COLORS: Record<MethodName, string> = {
  dirichlet: "#B06608", // chocolate
  parsimony: "#943F22",
  hss: "#09255D",
};

function methodTitle(m: MethodName) {
  return m === "dirichlet" ? "Dirichlet" : m === "parsimony" ? "Parsimony" : "HSS";
}

// Map ll_final values to visible pixel heights (handles negative/identical values)
function makeHeightScale(allLL: number[]) {
  const valid = allLL.filter((v) => Number.isFinite(v));
  if (!valid.length) {
    return (_v: number | null) => 80; // default
  }
  const min = Math.min(...valid);
  const max = Math.max(...valid);
  const span = (max - min) || 1;
  const minPx = 50;
  const maxPx = 180;
  return (v: number | null) => {
    if (!Number.isFinite(Number(v))) return minPx;
    const t = (Number(v) - min) / span;
    return Math.round(minPx + t * (maxPx - minPx));
  };
}

/** Build CSS repeating gradient for vertical stripes (equal white & color area). */
function stripeBackground(color: string, stripes: number) {
  const n = Math.max(0, stripes | 0);
  if (n === 0) return color;
  // One repeat spans 1/n of width; white occupies half, color the other half.
  // e.g., repeating-linear-gradient(90deg, white 0, white 50%/n, color 50%/n, color 100%/n)
  const whiteEnd = `calc(50% / ${n})`;
  const colorEnd = `calc(100% / ${n})`;
  return `repeating-linear-gradient(90deg, white 0, white ${whiteEnd}, ${color} ${whiteEnd}, ${color} ${colorEnd})`;
}

export default function CakeChart({ data, patternCounts, rep, root }: Props) {
  // bottom = largest cross-section = fine (layer 2), then medium (1), top = coarse (0)
  const order: (0 | 1 | 2)[] = [2, 1, 0];

  // precompute heights across all methods/layers for consistent scaling
  const allLL = useMemo(
    () =>
      (["parsimony", "dirichlet", "hss"] as MethodName[])
        .flatMap((m) => (data[m] || []))
        .map((L) => Number(L.ll_final))
        .filter((x) => Number.isFinite(x)),
    [data]
  );
  const height = useMemo(() => makeHeightScale(allLL), [allLL]);

  // cross-section width percent by layer from cumulative pattern weights
  const widthPctByLayer: Record<0 | 1 | 2, number> = {
    0: Math.max(0, Math.min(100, patternCounts.coarse)),
    1: Math.max(0, Math.min(100, patternCounts.medium)),
    2: Math.max(0, Math.min(100, patternCounts.fine)),
  };

  const methods: MethodName[] = ["parsimony", "dirichlet", "hss"];

  return (
    <div className="w-full">
      <div className="mb-3 text-sm text-white/90">
        {rep != null && <span className="mr-3">Rep: <b>{rep}</b></span>}
        {root && <span>Root: <b>{root}</b></span>}
      </div>

      <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
        {methods.map((method) => {
          const color = COLORS[method];
          const layers = (data[method] || []);
          // Index by layer id for quick lookup
          const byId: Partial<Record<0 | 1 | 2, OneLayer>> = {};
          for (const L of layers) byId[L.layer] = L;

          return (
            <div key={method} className="flex flex-col items-center">
              <div className="text-center mb-2 font-semibold tracking-wide">
                {methodTitle(method)}
              </div>

              {/* Cake stack */}
              <div className="relative w-full max-w-xs">
                {order.map((layerId, idx) => {
                  const L = byId[layerId] || { layer: layerId, iter: null, ll_final: null };
                  const hpx = height(L.ll_final);
                  const wPct = widthPctByLayer[layerId]; // 0..100
                  const bg = stripeBackground(color, Number(L.iter) || 0);

                  // center each tier; add small offset so upper tiers look stacked
                  return (
                    <div
                      key={layerId}
                      className="mx-auto rounded-md shadow-lg border border-white/40"
                      style={{
                        height: `${hpx}px`,
                        width: `${wPct}%`,
                        backgroundImage: bg,
                        backgroundColor: color,
                        backgroundSize: "auto 100%",
                        // tiny perspective illusion
                        transform: `translateY(${idx === 0 ? 0 : -6 * idx}px)`,
                      }}
                      title={`Layer ${layerId} · iter=${L.iter ?? 0} · ll_final=${L.ll_final ?? "-"}`}
                    />
                  );
                })}
              </div>
            </div>
          );
        })}
      </div>
    </div>
  );
}
