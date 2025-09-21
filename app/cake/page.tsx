// app/cake/page.tsx
"use client";

import React, { useEffect, useMemo, useRef, useState } from "react";
import { Canvas } from "@react-three/fiber";
import { OrthographicCamera, Edges, OrbitControls, Text } from "@react-three/drei";
import * as THREE from "three";

/* ===================== LAYER COLORS (EDIT IN CODE) ===================== */
/** Per-method colors for each layer: fine (bottom), medium (middle), coarse (top) */
const LAYER_COLOR = {
  parsimony: {
    fine:   "#EC5800", // persimmon
    medium: "#0D5EAF", // Greek blue
    coarse: "#FFFFFF", // white
  },
  hss: {
    fine:   "#7c3aed",
    medium: "#6d28d9",
    coarse: "#5b21b6",
  },
  dirichlet: {
    fine:   "#ef4444",
    medium: "#dc2626",
    coarse: "#b91c1c",
  },
} as const;
/* ====================================================================== */

/* ===================== DECORATIONS (EDIT IN CODE) ===================== */
/**
 * Supported kinds:
 *  - "diamonds": concentric diamond outlines
 *
 * placement:
 *  - "top" | "front" | "right" | "left"
 */
type DecorDiamond = {
  kind: "diamonds";
  placement?: "top" | "front" | "right" | "left";
  count?: number;
  color?: string;
  strokeRel?: number;
  insetRel?: number;
};

type DecorItem = DecorDiamond;
type LayerKey = "fine" | "medium" | "coarse";

/** Example: concentric diamonds on the RIGHT face of Parsimony's **middle** layer */
const DECOR_CONFIG: {
  parsimony?: Partial<Record<LayerKey, DecorItem[]>>;
  hss?: Partial<Record<LayerKey, DecorItem[]>>;
  dirichlet?: Partial<Record<LayerKey, DecorItem[]>>;
} = {
  parsimony: {
    medium: [
      { kind: "diamonds", placement: "right", count: 6, color: "#ffffff", strokeRel: 0.02, insetRel: 0.1 },
    ],
  },
  hss: {},
  dirichlet: {},
};
/* ====================================================================== */

/* ===================== HEIGHT SCALING (LOG OF RATIOS) ===================== */
const TARGET_MAX_CAKE_HEIGHT = 1.8;
const FALLBACK_SF = 0.004;
const LOG_EPS = 1e-9;
const MIN_LAYER_HEIGHT = 0.0001;

// ln(-x) with guards (x is negative LL)
const lnNeg = (x: number | null) => Math.log(Math.max(LOG_EPS, -Number(x ?? 0) || LOG_EPS));
// ln(|a|/|b|) with guards
const lnRatio = (a: number | null, b: number | null) => {
  const av = Math.max(LOG_EPS, Math.abs(Number(a ?? 0)));
  const bv = Math.max(LOG_EPS, Math.abs(Number(b ?? 0)));
  return Math.log(av / bv);
};

// bottom = ln(f/m), middle = ln(m/c), top = ln(-c)   (Total = ln(-f))
function heightsFromLL_logRatios(
  ll: { fine: number | null; medium: number | null; coarse: number | null },
  sf: number
): [number, number, number] {
  const b = Math.max(MIN_LAYER_HEIGHT, lnRatio(ll.fine, ll.medium) * sf);
  const m = Math.max(MIN_LAYER_HEIGHT, lnRatio(ll.medium, ll.coarse) * sf);
  const t = Math.max(MIN_LAYER_HEIGHT, lnNeg(ll.coarse) * sf);
  return [b, m, t];
}

type MethodKey = "parsimony" | "dirichlet" | "hss";
type CakeLL = { coarse: number | null; medium: number | null; fine: number | null };

function autoScaleFromTallest(llMap: Partial<Record<MethodKey, CakeLL>>): number {
  const totals: number[] = [];
  (["parsimony", "hss", "dirichlet"] as MethodKey[]).forEach((k) => {
    const ll = llMap[k];
    if (!ll) return;
    totals.push(lnNeg(ll.fine));
  });
  const maxTotal = Math.max(0, ...totals);
  return maxTotal > 0 ? TARGET_MAX_CAKE_HEIGHT / maxTotal : FALLBACK_SF;
}

/* ===================== PATTERN AREA SCALING ===================== */
/**
 * Convert site-pattern weights to world **area**.
 * side^2 (world units^2) = PATTERN_AREA_SF * weight
 */
const PATTERN_AREA_SF = 0.0036;  // <— EDIT to change layer cross-section size
const MIN_SIDE = 0.06;

function sideFromWeight(w: number | null | undefined): number {
  const area = Math.max(0, (Number(w ?? 0) || 0) * PATTERN_AREA_SF);
  return Math.max(MIN_SIDE, Math.sqrt(area));
}

/* ===================== ORIENTATION / LAYOUT ===================== */
const deg2rad = (d: number) => (d * Math.PI) / 180;

const cake_deg = -29;      // cake yaw (deg)
const placard_deg = -21;   // placard yaw (deg)
const CAKE_YAW_DEG = { parsimony: cake_deg, hss: cake_deg, dirichlet: cake_deg } as const;
const PLACARD_YAW_DEG = { parsimony: placard_deg, hss: placard_deg, dirichlet: placard_deg } as const;

const ROLL = Math.PI / 6;
const STAGE_SHIFT_Y = -2.1;
const X_BASE_SPREAD = 2.2;
const X_EXTRA_SEP = 0.9;
const X_SPREAD = X_BASE_SPREAD + X_EXTRA_SEP;
const BASELINE_Y = 0;

const PLACARD_CLEARANCE = 0.14;
const PLACARD_Y_DROP = 0.0;
const PLACARD_YAW_TWEAK = -Math.PI / 20;
const SHOW_HELPERS = false;

/* ===================== FETCH TYPES ===================== */
type PatternCounts = { coarse?: number; medium?: number; fine?: number };
type CakeResp = {
  ok?: boolean;
  job_id?: string;
  rep?: number;
  root?: string;
  ll: Partial<Record<MethodKey, CakeLL>>;
  pattern_counts?: PatternCounts;
};

type RepsResponse = { ok?: boolean; job_id?: string; reps: number[] };
type RootsResponse = { ok?: boolean; job_id?: string; rep?: number; roots: string[] };

/* ===================== Camera ===================== */
function IsoCamera({
  pos = [6, 7, 6],
  roll = ROLL,
  zoom = 200,
}: { pos?: [number, number, number]; roll?: number; zoom?: number }) {
  const ref = useRef<THREE.OrthographicCamera | null>(null);
  useEffect(() => {
    const cam = ref.current;
    if (!cam) return;
    cam.position.set(...(pos as [number, number, number]));
    cam.zoom = zoom;
    cam.updateProjectionMatrix();
    cam.lookAt(0, 0, 0);
    cam.rotation.z = roll;
  }, [pos, roll, zoom]);
  return <OrthographicCamera ref={ref} makeDefault near={0.01} far={1000} />;
}

/* ===================== Rounded-rectangle EXTRUDE (Option B) ===================== */
const SEAM_EPS = 0.0005;

/** Build a rounded-rectangle Shape centered at (0,0) in the XY plane. */
function roundedRectShape(width: number, height: number, radius: number): THREE.Shape {
  const w = Math.max(0, width);
  const h = Math.max(0, height);
  const r = Math.max(0, Math.min(radius, Math.min(w, h) / 2));

  const x = -w / 2;
  const y = -h / 2;

  const s = new THREE.Shape();
  s.moveTo(x + r, y);
  s.lineTo(x + w - r, y);
  s.quadraticCurveTo(x + w, y, x + w, y + r);
  s.lineTo(x + w, y + h - r);
  s.quadraticCurveTo(x + w, y + h, x + w - r, y + h);
  s.lineTo(x + r, y + h);
  s.quadraticCurveTo(x, y + h, x, y + h - r);
  s.lineTo(x, y + r);
  s.quadraticCurveTo(x, y, x + r, y);
  return s;
}

/* ===================== Decorations (textures & overlays) ===================== */
const DECOR_EPS = 0.001; // lift overlay slightly to avoid z-fighting

/** Square diamonds texture */
function makeDiamondTextureSquare({
  px = 1024,
  count = 5,
  color = "#ffffff",
  strokePx = 16,
  insetPx = 100,
}: {
  px?: number;
  count?: number;
  color?: string;
  strokePx?: number;
  insetPx?: number;
}): HTMLCanvasElement {
  const cvs = document.createElement("canvas");
  cvs.width = cvs.height = px;
  const ctx = cvs.getContext("2d")!;
  ctx.clearRect(0, 0, px, px);
  ctx.strokeStyle = color;
  ctx.lineWidth = Math.max(1, strokePx);
  ctx.lineJoin = "miter";
  ctx.lineCap = "butt";

  const cx = px / 2;
  const cy = px / 2;
  const outer = px / 2 - insetPx;
  const inner = Math.max(outer * 0.08, strokePx * 2);
  const steps = Math.max(1, count);

  for (let i = 0; i < steps; i++) {
    const t = steps === 1 ? 1 : i / (steps - 1);
    const r = inner + (outer - inner) * t; // half-diagonal length
    ctx.beginPath();
    ctx.moveTo(cx, cy - r);
    ctx.lineTo(cx + r, cy);
    ctx.lineTo(cx, cy + r);
    ctx.lineTo(cx - r, cy);
    ctx.closePath();
    ctx.stroke();
  }
  return cvs;
}

/** Diamonds texture for rectangular area (we keep it square; plane scales it) */
function makeDiamondTextureRect({
  width,
  height,
  count,
  color,
  strokeRel,
  insetRel,
  px = 1024,
}: {
  width: number; height: number; count?: number; color?: string; strokeRel?: number; insetRel?: number; px?: number;
}): THREE.CanvasTexture {
  const strokePx = Math.max(1, Math.round((strokeRel ?? 0.02) * px));
  const insetPx = Math.max(1, Math.round((insetRel ?? 0.1) * (px / 2)));

  const base = makeDiamondTextureSquare({
    px,
    count: count ?? 5,
    color: color ?? "#ffffff",
    strokePx,
    insetPx,
  });

  const texOut = new THREE.CanvasTexture(base);
  texOut.anisotropy = 8;
  texOut.wrapS = texOut.wrapT = THREE.ClampToEdgeWrapping;
  texOut.needsUpdate = true;
  return texOut;
}

/** Draw overlay on a chosen face */
function DecorOverlay({
  placement,
  side,
  height,
  centerY,
  items,
}: {
  placement: "top" | "front" | "right" | "left";
  side: number;
  height: number;     // used for vertical faces
  centerY: number;    // used for vertical faces
  items: DecorItem[] | undefined;
}) {
  const tex = useMemo(() => {
    if (!items || items.length === 0) return null;
    let canvasTex: THREE.CanvasTexture | null = null;
    for (const it of items) {
      if (it.kind === "diamonds") {
        const t = makeDiamondTextureRect({
          width: side,
          height: placement === "top" ? side : height,
          count: it.count,
          color: it.color,
          strokeRel: it.strokeRel,
          insetRel: it.insetRel,
        });
        canvasTex = t;
      }
    }
    return canvasTex;
  }, [items, placement, side, height]);

  if (!items || items.length === 0 || !tex) return null;

  if (placement === "top") {
    const topY = centerY + height / 2;
    return (
      <mesh position={[0, topY + DECOR_EPS, 0]} rotation={[-Math.PI / 2, 0, 0]}>
        <planeGeometry args={[side, side]} />
        <meshBasicMaterial map={tex} transparent opacity={1} side={THREE.DoubleSide} />
      </mesh>
    );
  }

  if (placement === "front") {
    return (
      <mesh position={[0, centerY, side / 2 + DECOR_EPS]}>
        <planeGeometry args={[side, height]} />
        <meshBasicMaterial map={tex} transparent opacity={1} side={THREE.DoubleSide} />
      </mesh>
    );
  }

  if (placement === "right") {
    return (
      <mesh position={[side / 2 + DECOR_EPS, centerY, 0]} rotation={[0, Math.PI / 2, 0]}>
        <planeGeometry args={[side, height]} />
        <meshBasicMaterial map={tex} transparent opacity={1} side={THREE.DoubleSide} />
      </mesh>
    );
  }

  // placement === "left"
  return (
    <mesh position={[-side / 2 - DECOR_EPS, centerY, 0]} rotation={[0, -Math.PI / 2, 0]}>
      <planeGeometry args={[side, height]} />
      <meshBasicMaterial map={tex} transparent opacity={1} side={THREE.DoubleSide} />
    </mesh>
  );
}

/**
 * Blunted cake layer using ExtrudeGeometry + optional decorations
 */
function BluntedLayer({
  side,
  height,
  color,
  centerY,
  yaw,
  topDecor,
  frontDecor,
  rightDecor,
  leftDecor,
  showEdges = true,
}: {
  side: number;
  height: number;
  color: string;
  centerY: number;
  yaw: number;
  topDecor?: DecorItem[] | undefined;
  frontDecor?: DecorItem[] | undefined;
  rightDecor?: DecorItem[] | undefined;
  leftDecor?: DecorItem[] | undefined;
  showEdges?: boolean;
}) {
  const cornerRadius = Math.min(side, height) * 0.06;

  const geo = useMemo(() => {
    const shape = roundedRectShape(side, side, cornerRadius);
    const g = new THREE.ExtrudeGeometry(shape, {
      depth: height,
      bevelEnabled: false,
      steps: 1,
    });
    g.rotateX(-Math.PI / 2);                 // make height along +Y
    g.translate(0, centerY - height / 2, 0); // bottom sits at baseline when centerY = baseline + h/2
    return g;
  }, [side, height, centerY, cornerRadius]);

  return (
    <group rotation={[0, yaw, 0]}>
      {/* Body */}
      <mesh geometry={geo}>
        <meshBasicMaterial color={color} />
        {showEdges && <Edges threshold={1} color="white" />}
      </mesh>

      {/* Decorations */}
      {topDecor && topDecor.length > 0 && (
        <DecorOverlay placement="top" side={side} height={height} centerY={centerY} items={topDecor} />
      )}
      {frontDecor && frontDecor.length > 0 && (
        <DecorOverlay placement="front" side={side} height={height} centerY={centerY} items={frontDecor} />
      )}
      {rightDecor && rightDecor.length > 0 && (
        <DecorOverlay placement="right" side={side} height={height} centerY={centerY} items={rightDecor} />
      )}
      {leftDecor && leftDecor.length > 0 && (
        <DecorOverlay placement="left" side={side} height={height} centerY={centerY} items={leftDecor} />
      )}
    </group>
  );
}

/** Helper type so component props accept readonly tuples */
type TriDecor = Readonly<[DecorItem[] | undefined, DecorItem[] | undefined, DecorItem[] | undefined]>;

/** Stacked cake using blunted layers (extruded rounded rectangles). */
function StackedCake({
  sides,                // [bottom (fine), middle (medium), top (coarse)]
  baselineY = BASELINE_Y,
  yaw,
  heights,             // [fineH, mediumH, coarseH]
  colors,              // [fineColor, mediumColor, coarseColor]
  decorsTop,           // triple of top decors
  decorsFront,         // triple of front decors
  decorsRight,         // triple of right decors
  decorsLeft,          // triple of left decors
}: {
  sides: [number, number, number];
  baselineY?: number;
  yaw: number;
  heights: [number, number, number];
  colors: [string, string, string];
  decorsTop?: TriDecor;
  decorsFront?: TriDecor;
  decorsRight?: TriDecor;
  decorsLeft?: TriDecor;
}) {
  const [sFine, sMed, sCoa] = sides;
  const [hFine, hMed, hCoa] = heights;
  const [cFine, cMed, cCoa] = colors;

  const [dTopFine, dTopMed, dTopCoa] = decorsTop ?? [undefined, undefined, undefined];
  const [dFrontFine, dFrontMed, dFrontCoa] = decorsFront ?? [undefined, undefined, undefined];
  const [dRightFine, dRightMed, dRightCoa] = decorsRight ?? [undefined, undefined, undefined];
  const [dLeftFine, dLeftMed, dLeftCoa] = decorsLeft ?? [undefined, undefined, undefined];

  const yFine = baselineY + hFine / 2;
  const yMed  = baselineY + hFine + SEAM_EPS + hMed / 2;
  const yCoa  = baselineY + hFine + SEAM_EPS + hMed + SEAM_EPS + hCoa / 2;

  return (
    <group>
      <BluntedLayer side={sFine} height={hFine} color={cFine} centerY={yFine} yaw={yaw}
        topDecor={dTopFine} frontDecor={dFrontFine} rightDecor={dRightFine} leftDecor={dLeftFine} />
      <BluntedLayer side={sMed}  height={hMed}  color={cMed}  centerY={yMed}  yaw={yaw}
        topDecor={dTopMed} frontDecor={dFrontMed} rightDecor={dRightMed} leftDecor={dLeftMed} />
      <BluntedLayer side={sCoa}  height={hCoa}  color={cCoa}  centerY={yCoa}  yaw={yaw}
        topDecor={dTopCoa} frontDecor={dFrontCoa} rightDecor={dRightCoa} leftDecor={dLeftCoa} />
    </group>
  );
}

function TentPlacard({
  label, baselineY = BASELINE_Y, yawForOffset, yaw, yNudge = 0,
  width = 0.7, height = 0.252, angleDeg = 26, textSize = 0.084, adjacentCubeSide,
}: {
  label: string; baselineY?: number; yawForOffset: number; yaw: number; yNudge?: number;
  width?: number; height?: number; angleDeg?: number; textSize?: number; adjacentCubeSide: number;
}) {
  const theta = (angleDeg * Math.PI) / 180;
  const hingeY = baselineY + height * Math.cos(theta) - PLACARD_Y_DROP + yNudge;
  const rx = Math.cos(yawForOffset), rz = -Math.sin(yawForOffset);
  const pos: [number, number, number] = [
    rx * (adjacentCubeSide / 2 + PLACARD_CLEARANCE),
    hingeY,
    rz * (adjacentCubeSide / 2 + PLACARD_CLEARANCE),
  ];
  const yRot = yaw + Math.PI / 2 + PLACARD_YAW_TWEAK;
  const matWhite = new THREE.MeshBasicMaterial({ color: "#ffffff", side: THREE.DoubleSide });

  return (
    <group position={pos} rotation={[0, yRot, 0]}>
      <group rotation={[-theta, 0, 0]}>
        <mesh position={[0, -height / 2, 0]} material={matWhite}>
          <planeGeometry args={[width, height]} />
        </mesh>
        <Text position={[0, -height / 2 + 0.045, 0.003]} color="black" fontSize={textSize} anchorX="center" anchorY="top" maxWidth={width * 0.92}>
          {label}
        </Text>
      </group>
      <group rotation={[theta, 0, 0]}>
        <mesh position={[0, -height / 2, 0]} material={matWhite}>
          <planeGeometry args={[width, height]} />
        </mesh>
      </group>
    </group>
  );
}

/* ===================== Page ===================== */
export default function CakesPage() {
  // selectors
  const [job, setJob] = useState<string>("");
  const [reps, setReps] = useState<number[]>([]);
  const [rep, setRep] = useState<number>(4);        // default
  const [roots, setRoots] = useState<string[]>([]);
  const [root, setRoot] = useState<string>("h_6"); // default

  // data: log-likelihoods + pattern counts
  const [llMap, setLlMap] = useState<Partial<Record<MethodKey, CakeLL>>>({});
  const [patternCounts, setPatternCounts] = useState<PatternCounts | null>(null);

  /* ---------- read job from URL or localStorage ---------- */
  useEffect(() => {
    try {
      const u = new URL(window.location.href);
      const fromUrl = u.searchParams.get("job") || u.searchParams.get("job_id") || "";
      const fromLs = localStorage.getItem("emtr:selectedJobId") || "";
      const pick = (fromUrl || fromLs || "").trim();
      if (pick) setJob(pick);
    } catch { /* ignore */ }
  }, []);

  /* ---------- reps ---------- */
  useEffect(() => {
    (async () => {
      if (!job) return;
      setReps([]);
      setRoots([]);
      setRoot("h_6");
      try {
        const r = await fetch(`/api/results/${encodeURIComponent(job)}/reps`, { cache: "no-store" });
        const j: RepsResponse = await r.json();
        const list = (Array.isArray(j.reps) ? j.reps : []).filter((n) => Number.isFinite(n));
        setReps(list);
        if (!list.includes(rep) && list.length) setRep(list[0]);
      } catch { /* ignore */ }
    })();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [job]);

  /* ---------- roots (depends on rep) ---------- */
  useEffect(() => {
    (async () => {
      if (!job || rep == null) return;
      setRoots([]);
      setRoot("h_6");
      try {
        const r = await fetch(`/api/results/${encodeURIComponent(job)}/reps/roots?rep=${rep}`, { cache: "no-store" });
        const j: RootsResponse = await r.json();
        const list = Array.isArray(j.roots) ? j.roots : [];
        setRoots(list);
        if (list.length) setRoot((prev) => (prev && list.includes(prev) ? prev : list.includes("h_6") ? "h_6" : list[0]));
      } catch { /* ignore */ }
    })();
  }, [job, rep]);

  /* ---------- cake LL & pattern counts ---------- */
  useEffect(() => {
    (async () => {
      if (!job || rep == null || !root) return;

      try {
        // 1) LLs + pattern weights
        const r1 = await fetch(`/api/results/${encodeURIComponent(job)}/cake?rep=${rep}&root=${encodeURIComponent(root)}`, { cache: "no-store" });
        const j1: CakeResp = await r1.json();
        if (j1?.ll) setLlMap(j1.ll);
        setPatternCounts(j1?.pattern_counts ?? null);
      } catch { /* ignore */ }
    })();
  }, [job, rep, root]);

  // Compute scale factor (auto) and derived sides/heights per method
  const sf = useMemo(() => autoScaleFromTallest(llMap), [llMap]);

  // pattern weights → sides (fine bottom, medium middle, coarse top)
  const sidesByMethod = useMemo(() => {
    const pc = patternCounts || { fine: 100, medium: 60, coarse: 30 };
    const bottomSide = sideFromWeight(pc.fine);
    const middleSide = sideFromWeight(pc.medium);
    const topSide = sideFromWeight(pc.coarse);
    return {
      parsimony: [bottomSide, middleSide, topSide] as [number, number, number],
      hss: [bottomSide, middleSide, topSide] as [number, number, number],
      dirichlet: [bottomSide, middleSide, topSide] as [number, number, number],
    };
  }, [patternCounts]);

  // heights per method (log-ratio definition)
  const heightsByMethod = useMemo(() => {
    const out: Record<MethodKey, [number, number, number]> = {
      parsimony: [0.1, 0.1, 0.1],
      hss: [0.1, 0.1, 0.1],
      dirichlet: [0.1, 0.1, 0.1],
    };
    (["parsimony", "hss", "dirichlet"] as MethodKey[]).forEach((k) => {
      const ll = llMap[k];
      if (!ll) return;
      out[k] = heightsFromLL_logRatios(ll, sf);
    });
    return out;
  }, [llMap, sf]);

  // colors per method → [fine, medium, coarse]
  const colorsByMethod = useMemo(() => {
    const toTuple = (m: keyof typeof LAYER_COLOR) =>
      [LAYER_COLOR[m].fine, LAYER_COLOR[m].medium, LAYER_COLOR[m].coarse] as [string, string, string];
    return {
      parsimony: toTuple("parsimony"),
      hss: toTuple("hss"),
      dirichlet: toTuple("dirichlet"),
    };
  }, []);

  // decorations per method, grouped by placement
  const pickDecor = (m: "parsimony" | "hss" | "dirichlet", layer: LayerKey, placement: DecorDiamond["placement"]) =>
    (DECOR_CONFIG[m] && DECOR_CONFIG[m]![layer]?.filter(d => (d.placement ?? "top") === placement)) || undefined;

  type TriDecor = Readonly<[DecorItem[] | undefined, DecorItem[] | undefined, DecorItem[] | undefined]>;

  const decorsTopByMethod: Record<MethodKey, TriDecor> = useMemo(() => ({
    parsimony: [pickDecor("parsimony", "fine", "top"), pickDecor("parsimony", "medium", "top"), pickDecor("parsimony", "coarse", "top")] as const,
    hss:       [pickDecor("hss", "fine", "top"),       pickDecor("hss", "medium", "top"),       pickDecor("hss", "coarse", "top")] as const,
    dirichlet: [pickDecor("dirichlet", "fine", "top"), pickDecor("dirichlet", "medium", "top"), pickDecor("dirichlet", "coarse", "top")] as const,
  }), []);

  const decorsFrontByMethod: Record<MethodKey, TriDecor> = useMemo(() => ({
    parsimony: [pickDecor("parsimony", "fine", "front"), pickDecor("parsimony", "medium", "front"), pickDecor("parsimony", "coarse", "front")] as const,
    hss:       [pickDecor("hss", "fine", "front"),       pickDecor("hss", "medium", "front"),       pickDecor("hss", "coarse", "front")] as const,
    dirichlet: [pickDecor("dirichlet", "fine", "front"), pickDecor("dirichlet", "medium", "front"), pickDecor("dirichlet", "coarse", "front")] as const,
  }), []);

  const decorsRightByMethod: Record<MethodKey, TriDecor> = useMemo(() => ({
    parsimony: [pickDecor("parsimony", "fine", "right"), pickDecor("parsimony", "medium", "right"), pickDecor("parsimony", "coarse", "right")] as const,
    hss:       [pickDecor("hss", "fine", "right"),       pickDecor("hss", "medium", "right"),       pickDecor("hss", "coarse", "right")] as const,
    dirichlet: [pickDecor("dirichlet", "fine", "right"), pickDecor("dirichlet", "medium", "right"), pickDecor("dirichlet", "coarse", "right")] as const,
  }), []);

  const decorsLeftByMethod: Record<MethodKey, TriDecor> = useMemo(() => ({
    parsimony: [pickDecor("parsimony", "fine", "left"), pickDecor("parsimony", "medium", "left"), pickDecor("parsimony", "coarse", "left")] as const,
    hss:       [pickDecor("hss", "fine", "left"),       pickDecor("hss", "medium", "left"),       pickDecor("hss", "coarse", "left")] as const,
    dirichlet: [pickDecor("dirichlet", "fine", "left"), pickDecor("dirichlet", "medium", "left"), pickDecor("dirichlet", "coarse", "left")] as const,
  }), []);

  // positions/orientation helpers
  const ux = Math.cos(ROLL);
  const uy = Math.sin(ROLL);

  // yaws
  const cakeYawLeft = deg2rad(CAKE_YAW_DEG.parsimony);
  const cakeYawCenter = deg2rad(CAKE_YAW_DEG.hss);
  const cakeYawRight = deg2rad(CAKE_YAW_DEG.dirichlet);
  const placYawLeft = deg2rad(PLACARD_YAW_DEG.parsimony);
  const placYawCenter = deg2rad(PLACARD_YAW_DEG.hss);
  const placYawRight = deg2rad(PLACARD_YAW_DEG.dirichlet);

  // baselines
  const basePar = BASELINE_Y;
  const baseHss = BASELINE_Y;
  const baseDir = BASELINE_Y;

  // pick adjacentCubeSide = bottom side for each method
  const sPar = sidesByMethod.parsimony[0];
  const sHss = sidesByMethod.hss[0];
  const sDir = sidesByMethod.dirichlet[0];

  return (
    <div className="relative min-h-screen">
      {/* Background */}
      <div
        className="absolute inset-0 -z-10"
        style={{
          backgroundImage: `url(/assets/evergreen.png)`,
          backgroundSize: "cover",
          backgroundPosition: "center",
          filter: "contrast(90%) brightness(115%) blur(0.3px)",
          opacity: 0.75,
        }}
      />

      {/* Controls */}
      <div className="relative z-10 max-w-5xl mx-auto p-4 flex flex-wrap gap-3 items-end text-sm text-black">
        <div>Job: <span className="font-mono">{job || "(none)"} </span></div>
        <label className="ml-4">
          Rep:&nbsp;
          <select
            className="border rounded px-2 py-1"
            value={rep}
            onChange={(e) => setRep(Number(e.target.value))}
            disabled={!reps.length}
          >
            {reps.map((r) => <option key={r} value={r}>{r}</option>)}
          </select>
        </label>
        <label>
          Root:&nbsp;
          <select
            className="border rounded px-2 py-1"
            value={root}
            onChange={(e) => setRoot(e.target.value)}
            disabled={!roots.length}
          >
            {roots.map((r) => <option key={r} value={r}>{r}</option>)}
          </select>
        </label>
        {!job && (
          <div className="text-yellow-700">
            Select a job on the Available Results page — or add ?job=... to the URL.
          </div>
        )}
      </div>

      <Canvas gl={{ antialias: true, alpha: true }} dpr={[1, 2]} style={{ position: "absolute", inset: 0 }}>
        <IsoCamera pos={[6, 7, 6]} roll={ROLL} zoom={200} />

        {/* Stage near bottom; all stacks share the same screen-aligned baseline */}
        <group position={[0, STAGE_SHIFT_Y, 0]}>
          {SHOW_HELPERS && (
            <>
              <gridHelper args={[12, 12, 0x888888, 0x444444]} position={[0, BASELINE_Y - 0.001, 0]} />
              <axesHelper args={[3]} />
            </>
          )}

          {/* LEFT: Parsimony */}
          <group position={[-X_SPREAD * ux, -X_SPREAD * uy, 0]}>
            <StackedCake
              sides={sidesByMethod.parsimony}
              heights={heightsByMethod.parsimony}
              colors={colorsByMethod.parsimony}
              decorsTop={decorsTopByMethod.parsimony}
              decorsFront={decorsFrontByMethod.parsimony}
              decorsRight={decorsRightByMethod.parsimony}
              decorsLeft={decorsLeftByMethod.parsimony}
              baselineY={basePar}
              yaw={cakeYawLeft}
            />
            <TentPlacard
              label="Parsimony"
              baselineY={basePar}
              adjacentCubeSide={sPar}
              yawForOffset={cakeYawLeft}
              yaw={placYawLeft}
            />
          </group>

          {/* CENTER: HSS */}
          <group position={[0, 0, 0]}>
            <StackedCake
              sides={sidesByMethod.hss}
              heights={heightsByMethod.hss}
              colors={colorsByMethod.hss}
              decorsTop={decorsTopByMethod.hss}
              decorsFront={decorsFrontByMethod.hss}
              decorsRight={decorsRightByMethod.hss}
              decorsLeft={decorsLeftByMethod.hss}
              baselineY={baseHss}
              yaw={cakeYawCenter}
            />
            <TentPlacard
              label="HSS"
              baselineY={baseHss}
              adjacentCubeSide={sHss}
              yawForOffset={cakeYawCenter}
              yaw={placYawCenter}
            />
          </group>

          {/* RIGHT: Dirichlet */}
          <group position={[X_SPREAD * ux, X_SPREAD * uy, 0]}>
            <StackedCake
              sides={sidesByMethod.dirichlet}
              heights={heightsByMethod.dirichlet}
              colors={colorsByMethod.dirichlet}
              decorsTop={decorsTopByMethod.dirichlet}
              decorsFront={decorsFrontByMethod.dirichlet}
              decorsRight={decorsRightByMethod.dirichlet}
              decorsLeft={decorsLeftByMethod.dirichlet}
              baselineY={baseDir}
              yaw={cakeYawRight}
            />
            <TentPlacard
              label="Dirichlet"
              baselineY={baseDir}
              adjacentCubeSide={sDir}
              yawForOffset={cakeYawRight}
              yaw={placYawRight}
            />
          </group>
        </group>

        <OrbitControls enableRotate={false} enablePan={false} enableZoom={false} />
      </Canvas>
    </div>
  );
}
