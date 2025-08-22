// /lib/exportName.ts
type ExportNameParams = {
  job?: string;                 // e.g. wasm-17558...
  root?: string;                // violin page
  thr?: number;
  reps?: number;
  maxIter?: number;
  D_pi?: number[];              // [a,b,c,d]
  D_M?: number[];               // [a,b,c,d]
  twistDeg?: number;            // spiral pose
  tag?: string;                 // optional extra tag you might want
};

function sanitize(s: string) {
  return String(s).replace(/[/\\?%*:|"<>]/g, "-");
}

export function buildExportPrefix(p: ExportNameParams): string {
  const parts: string[] = [];

  if (p.job) parts.push(sanitize(p.job));
  if (p.root) parts.push(`root${sanitize(p.root)}`);
  if (Number.isFinite(p.thr as number)) parts.push(`thr${p.thr}`);
  if (Number.isFinite(p.reps as number)) parts.push(`reps${p.reps}`);
  if (Number.isFinite(p.maxIter as number)) parts.push(`maxIter${p.maxIter}`);
  if (Array.isArray(p.D_pi) && p.D_pi.length === 4) parts.push(`Dpi-${p.D_pi.join("-")}`);
  if (Array.isArray(p.D_M)  && p.D_M.length  === 4) parts.push(`DM-${p.D_M.join("-")}`);
  if (Number.isFinite(p.twistDeg as number)) parts.push(`tw${Math.round(p.twistDeg as number)}deg`);
  if (p.tag) parts.push(sanitize(p.tag));

  // Always return *something* useful.
  return parts.length ? parts.join("__") : "export";
}
