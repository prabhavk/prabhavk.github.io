// public/wasm/worker.js

// ---------- Upload config ----------
const API_BASE = self.API_BASE || "";
let AUTH = "";
let API_SECRET = "";

function buildAuthHeaders() {
  const h = { "content-type": "application/json" };
  if (AUTH) h.authorization = AUTH;
  if (API_SECRET) h["x-emtr-secret"] = API_SECRET;
  return h;
}

// Allow overriding endpoint paths without rebuilding the worker
const ENDPOINTS = {
  rows:    (jobId) => `${API_BASE}/api/jobs/${encodeURIComponent(jobId)}/rows`,
  allinfo: () => `${API_BASE}${self.ALL_INFO_PATH || "/api/allinfo"}`,
  trees:   () => `${API_BASE}/api/trees`,
  aa:      () => `${API_BASE}/api/aa-info`,
};

// ---------- TOGGLES ----------
let DB_LOGS = typeof self.DB_LOGS === "boolean" ? self.DB_LOGS : false;
let UPSERT_UI_LOGS = typeof self.UPSERT_UI_LOGS === "boolean" ? self.UPSERT_UI_LOGS : false;
let UPLOAD_TREES = typeof self.UPLOAD_TREES === "boolean" ? self.UPLOAD_TREES : true;

// Helpers for logging
function dbLog(line) {
  if (!DB_LOGS) return;
  try {
    const ts = new Date().toISOString();
    console.debug(`[DB][${ts}][${jobIdForThisRun || "-"}] ${line}`);
  } catch {}
}
const upLog = (line) => { if (UPSERT_UI_LOGS) postMessage({ type: "log", line }); };

// batching/backoff
const MAX_BATCH = 200;
const FLUSH_INTERVAL_MS = 1000;
const RETRY_BASE_MS = 1500;
const RETRY_MAX_MS = 15000;

// Beacon payload budget (conservative; some browsers allow ~64KB)
const BEACON_MAX_BYTES = 60000;

// ---------- Upload state ----------
let jobIdForThisRun = "";
let rowBuffer = [];
let flushTimer = null;
let inflight = false;
let retryDelay = RETRY_BASE_MS;
let startedAtMs = 0;
let flushSeq = 0;

// rows endpoint circuit breaker
let rowsEndpointDisabled = false;
let rowsEndpointErrorReason = "";
let artifactSeq = 0;

// Track latest repetition number seen in streaming [ROW {...}]
let lastRepSeen = null;

// ---------------- [EM_AllInfo] coalescing state (per-rep) ----------------
let bestFlushTimer = null;
let bestEndpointMissing = false;
let lastBestObj = null;
let bestUploaded = false;

const BEST_UPLOAD_DEBOUNCE_MS = 500;

// Keep a pending bundle for the *current* rep only
let pendingRep = null;
const pendingBest = { parsimony: null, dirichlet: null, hss: null };

// ---------------- Trees state ----------------
let treesUploadedCount = 0;   // mark first as current
let epsilonParam = null;      // captured from params if provided
const seenTreeHashesByJob = new Map(); // jobId -> Set(hash) for per-job dedupe
const pendingTreeUploads = new Set();  // in-flight tree uploads to await at end

// ---------- generic helpers ----------
function approxBytesOfRow(r) { return 8 + JSON.stringify(r).length + 1; }
function approxBytesOfJson(obj) { return JSON.stringify(obj).length; }
function canUseBeacon() { try { return typeof navigator !== "undefined" && typeof navigator.sendBeacon === "function"; } catch { return false; } }
function tryBeacon(url, obj) { try { const blob = new Blob([JSON.stringify(obj)], { type: "application/json" }); return navigator.sendBeacon(url, blob); } catch { return false; } }
async function postJSON(url, obj, { preferBeacon = false, keepalive = false } = {}) {
  if (preferBeacon && canUseBeacon() && approxBytesOfJson(obj) <= BEACON_MAX_BYTES) {
    const ok = tryBeacon(url, obj); if (ok) return { ok: true, via: "beacon" };
  }
  const res = await fetch(url, {
    method: "POST",
    headers: buildAuthHeaders(),
    body: JSON.stringify(obj),
    keepalive,
    cache: "no-store",
  });
  if (!res.ok) {
    const text = await res.text().catch(() => "");
    const err = new Error(`HTTP ${res.status} ${res.statusText} :: ${text.slice(0, 400)}`);
    err.status = res.status; err.bodyText = text || ""; err.permanent = res.status === 404 || res.status === 410;
    throw err;
  }
  return { ok: true, via: "fetch" };
}
function classifyPermanentRowsError(err) {
  if (!err) return { permanent: false, reason: "" };
  const status = err.status | 0;
  const msg = (err.bodyText || err.message || "").toLowerCase();
  if (status === 404 || status === 410) return { permanent: true, reason: `endpoint returned ${status}` };
  const tableMissing =
    msg.includes("doesn't exist") || msg.includes("does not exist") ||
    msg.includes("no such table") || msg.includes("unknown table") ||
    (msg.includes("table") && msg.includes("not found"));
  if (status === 500 && tableMissing) return { permanent: true, reason: "target table is missing on the server" };
  return { permanent: false, reason: "" };
}
function stripAnsi(s) { return String(s).replace(/\x1b\[[0-9;]*m/g, ""); }
function ensureSemicolon(nwk) { const s = String(nwk).trim(); return /;\s*$/.test(s) ? s : s + ";"; }
function safeNumber(x) { const n = Number(x); return Number.isFinite(n) ? n : null; }
function hashString(str) { let h = 5381, i = str.length; while (i) h = (h * 33) ^ str.charCodeAt(--i); return (h >>> 0).toString(36); }
function sanitize(name) { if (!name) return null; const c = String(name).replace(/[^A-Za-z0-9._-]/g, "_"); return c.length ? c : null; }

// ---------- Job metadata upsert ----------
async function upsertJobMetadata(jobId, cfg) {
  const payload = { job_id: jobId, thr: cfg.thr, reps: cfg.reps, max_iter: cfg.maxIter, D_pi: cfg.D_pi, D_M: cfg.D_M, status: "started" };
  const size = approxBytesOfJson(payload);
  try {
    upLog(`↗ upsert /api/jobs (job=${jobId}, ~${size}B)`);
    const result = await postJSON(`${API_BASE}/api/jobs`, payload, { keepalive: true });
    dbLog(`job config saved via ${result.via}`);
    upLog(`✔ upsert /api/jobs ok (via ${result.via})`);
  } catch (err) {
    postMessage({ type: "log", line: `❌ upsert /api/jobs failed: ${String(err)}` });
  }
}

// ---------- Job status patch ----------
async function setJobStatus(jobId, payload) {
  try {
    const res = await fetch(`${API_BASE}/api/jobs/${encodeURIComponent(jobId)}`, {
      method: "PATCH",
      headers: buildAuthHeaders(),
      keepalive: true,
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const text = await res.text().catch(() => "");
      postMessage({ type: "log", line: `⚠️ job patch failed: HTTP ${res.status} :: ${text.slice(0,200)}` });
    } else {
      if (payload?.status) dbLog(`job status -> ${payload.status}`);
    }
  } catch (e) {
    postMessage({ type: "log", line: `⚠️ job patch error: ${String(e)}` });
  }
}

function summarizeBatch(batch) {
  const methods = Object.create(null);
  let minRep = Infinity, maxRep = -Infinity;
  let minIter = Infinity, maxIter = -Infinity;
  let approxBytes = 0;
  for (const r of batch) {
    const m = (r.method ?? "(none)").toString();
    methods[m] = (methods[m] || 0) + 1;
    if (Number.isFinite(r.rep)) { minRep = Math.min(minRep, r.rep); maxRep = Math.max(maxRep, r.rep); }
    if (Number.isFinite(r.iter)) { minIter = Math.min(minIter, r.iter); maxIter = Math.max(maxIter, r.iter); }
    approxBytes += approxBytesOfRow(r);
  }
  if (!Number.isFinite(minRep)) { minRep = null; maxRep = null; }
  if (!Number.isFinite(minIter)) { minIter = null; maxIter = null; }
  const methodStr = Object.keys(methods).length
    ? Object.entries(methods).map(([k,v]) => `${k}:${v}`).join(", ")
    : "—";
  return { methodStr, minRep, maxRep, minIter, maxIter, approxBytes };
}

function normalizeRow(raw) {
  const toNum = (v) => (v == null || v === "" ? NaN : Number(v));
  const toInt = (v) => (v == null || v === "" ? NaN : parseInt(v, 10));
  const ll_init_src   = raw.ll_init ?? raw.ll_initial ?? raw.logLikelihood_initial ?? raw.loglikelihood_initial;
  const ecd_first_src = raw.ecd_ll_first ?? raw.logLikelihood_ecd_first ?? raw.loglikelihood_ecd_first;
  const ecd_final_src = raw.ecd_ll_final ?? raw.logLikelihood_ecd_final ?? raw.loglikelihood_ecd_final;
  const ll_final_src  = raw.ll_final ?? raw.logLikelihood_final ?? raw.loglikelihood_final;
  const rep_src  = raw.rep  ?? raw.repetition ?? raw.repeat;
  const iter_src = raw.iter ?? raw.iteration;
  const root_src = raw.root ?? raw.start ?? raw.node;
  return {
    method: raw.method,
    root: String(root_src ?? ""),
    rep: Number.isInteger(raw.rep) ? raw.rep : toInt(rep_src),
    iter: Number.isInteger(raw.iter) ? raw.iter : toInt(iter_src),
    ll_init: toNum(ll_init_src),
    ecd_ll_first: toNum(ecd_first_src),
    ecd_ll_final: toNum(ecd_final_src),
    ll_final: toNum(ll_final_src),
  };
}

function scheduleFlush() {
  if (!flushTimer) flushTimer = setTimeout(() => void flushNow(jobIdForThisRun, { final: false }), FLUSH_INTERVAL_MS);
}
function queueRow(jobId, row) {
  const r = normalizeRow(row);
  if (!r) return;
  if (Number.isFinite(r.rep)) lastRepSeen = r.rep;
  rowBuffer.push(r);
  if (rowBuffer.length >= MAX_BATCH) {
    if (!inflight) void flushNow(jobId, { final: false });
    else scheduleFlush();
  } else {
    scheduleFlush();
  }
}
function takeByteCappedBatch(maxRows, maxBytes) {
  const batch = [];
  let bytes = 0;
  while (batch.length < maxRows && rowBuffer.length) {
    const next = rowBuffer[0];
    const cost = approxBytesOfRow(next);
    if (batch.length > 0 && bytes + cost > maxBytes) break;
    batch.push(rowBuffer.shift());
    bytes += cost;
  }
  return batch;
}
function encodeCsv(rows) {
  const header = ["job_id","method","root","rep","iter","ll_init","ecd_ll_first","ecd_ll_final","ll_final"];
  const toCell = (v) => {
    if (v == null || Number.isNaN(v)) return "";
    const s = String(v);
    return /[",\n]/.test(s) ? `"${s.replace(/"/g,'""')}"` : s;
  };
  const lines = [header.join(",")];
  for (const r of rows) {
    lines.push([
      jobIdForThisRun,
      r.method ?? "",
      r.root ?? "",
      Number.isFinite(r.rep) ? r.rep : "",
      Number.isFinite(r.iter) ? r.iter : "",
      Number.isFinite(r.ll_init) ? r.ll_init : "",
      Number.isFinite(r.ecd_ll_first) ? r.ecd_ll_first : "",
      Number.isFinite(r.ecd_ll_final) ? r.ecd_ll_final : "",
      Number.isFinite(r.ll_final) ? r.ll_final : "",
    ].map(toCell).join(","));
  }
  return new TextEncoder().encode(lines.join("\n"));
}
function emitRowsAsArtifact(rows, reason) {
  try {
    const name = `llchange_fallback_${jobIdForThisRun}_${++artifactSeq}.csv`;
    const bytes = encodeCsv(rows);
    postMessage({ type: "log", line: `📦 Rows saved as artifact "${name}" (${bytes.length}B) due to: ${reason}` });
    postMessage({ type: "artifact", name, bytes }, [bytes.buffer]);
  } catch (e) {
    postMessage({ type: "log", line: `⚠️ failed to emit fallback artifact: ${String(e)}` });
  }
}

/* ---------------- Trees uploading helpers ---------------- */

function estimateCountsFromNewick(nw) {
  const n_edges = (nw.match(/:/g) || []).length || null;
  const leafLike = nw.match(/([^(),:]+)(?=:[^(),;]*)/g) || [];
  const n_leaves = leafLike.length || null;
  return { n_leaves, n_edges };
}

// Robustly parse [FJ_TREE] {json or raw} / [TREE] ... / NEWICK:
function maybeParseTreeLine(lineStr, jobId) {
  let line = stripAnsi(String(lineStr)).trim();
  if (!line) return null;

  const iFJ = line.indexOf("[FJ_TREE]");
  const iT  = iFJ === -1 ? line.indexOf("[TREE]") : -1;
  const tagAt = iFJ !== -1 ? iFJ : iT;
  if (tagAt === -1) {
    if (line.startsWith("NEWICK:")) {
      const nwk = line.slice("NEWICK:".length).trim();
      if (nwk) {
        const meta = {}; const jid = String(jobId ?? "").trim();
        if (jid) meta.job_id = jid;
        return { newick: ensureSemicolon(nwk), meta };
      }
    }
    return null;
  }

  line = line.slice(tagAt);
  const m = line.match(/^\[(?:FJ_)?TREE\]\s*(.+)$/);
  if (!m) return null;
  const payload = m[1].trim();
  if (!payload) return null;

  if (payload[0] === "{") {
    try {
      const obj = JSON.parse(payload);
      const newick =
        obj.newick != null ? String(obj.newick) :
        obj.tree   != null ? String(obj.tree)   :
        obj.nwk    != null ? String(obj.nwk)    : "";
      if (!newick) return null;

      const meta = {};
      const jid = String(obj.job_id ?? "").trim() || String(jobId ?? "").trim();
      if (jid) meta.job_id = jid;

      if (obj.label != null)     meta.label = String(obj.label);
      if (obj.method != null)    meta.method = String(obj.method);
      if (obj.root != null)      meta.root = String(obj.root);

      const rep  = Number(obj.rep);                if (Number.isFinite(rep))  meta.rep  = rep;
      const iter = Number(obj.iter);               if (Number.isFinite(iter)) meta.iter = iter;
      const eps  = Number(obj.epsilon ?? obj.eps); if (Number.isFinite(eps))  meta.epsilon = eps;

      return { newick: ensureSemicolon(newick), meta };
    } catch {}
  }

  const meta = {}; const jid2 = String(jobId ?? "").trim();
  if (jid2) meta.job_id = jid2;
  return { newick: ensureSemicolon(payload), meta };
}

async function postTree(jobId, newick, meta = {}) {
  if (!UPLOAD_TREES) return false;
  if (!jobId) { postMessage({ type: "log", line: "⚠️ postTree: missing jobId" }); return false; }
  if (!newick || typeof newick !== "string") return false;

  const nwk = ensureSemicolon(stripAnsi(newick));
  const jobSeen = seenTreeHashesByJob.get(jobId) || new Set();
  const dedupeKey = `${nwk}|${meta.label ?? ""}|${meta.rep ?? ""}|${meta.iter ?? ""}`;
  const h = hashString(dedupeKey);
  if (jobSeen.has(h)) { dbLog(`postTree: dedupe skip (hash=${h})`); return false; }

  let n_leaves = null, n_edges = null;
  try { const est = estimateCountsFromNewick(nwk); n_leaves = est.n_leaves ?? null; n_edges = est.n_edges ?? null; } catch {}

  const payload = {
    job_id: jobId,
    source: meta.source ?? "familyjoining",
    method: meta.method ?? "na",
    label:  meta.label  ?? "final",
    format: "newick",
    tree:   nwk,
    rep:    safeNumber(meta.rep),
    iter:   safeNumber(meta.iter),
    epsilon: safeNumber(meta.epsilon ?? (typeof epsilonParam === "number" ? epsilonParam : null)),
    n_leaves,
    n_edges,
    root: meta.root ?? null,
    is_current: meta.is_current ?? (treesUploadedCount === 0 ? 1 : 0),
  };

  try {
    upLog(`↗ upsert /api/trees (len=${nwk.length})`);
    const res = await postJSON(ENDPOINTS.trees(), payload, { keepalive: true });
    treesUploadedCount++;
    jobSeen.add(h);
    seenTreeHashesByJob.set(jobId, jobSeen);
    dbLog(`tree saved via ${res.via ?? "db"} (current=${payload.is_current})`);
    upLog(`✔ upsert /api/trees ok ${res.via ? `(via ${res.via})` : ""}`);
    return true;
  } catch (e) {
    postMessage({ type: "log", line: `⚠️ tree upload failed: ${String(e)}` });
    return false;
  }
}
async function awaitTreeUploads() {
  try {
    const all = Array.from(pendingTreeUploads);
    if (all.length) await Promise.allSettled(all);
  } catch {}
}

/* ---------------- rows flush ---------------- */

async function flushNow(jobId, { final }) {
  if (inflight) return false;
  if (flushTimer) { clearTimeout(flushTimer); flushTimer = null; }
  if (!rowBuffer.length) return false;

  if (rowsEndpointDisabled) {
    const batch = rowBuffer.splice(0, Math.min(rowBuffer.length, MAX_BATCH));
    const id = ++flushSeq;
    emitRowsAsArtifact(batch, rowsEndpointErrorReason || "rows endpoint disabled");
    postMessage({ type: "log", line: `[sent #${id}] skipped network (endpoint disabled) · remaining buffer: ${rowBuffer.length}` });
    if (rowBuffer.length) setTimeout(() => void flushNow(jobId, { final }), 0);
    return true;
  }

  const batch = final && canUseBeacon()
    ? takeByteCappedBatch(MAX_BATCH, BEACON_MAX_BYTES - 512)
    : rowBuffer.splice(0, Math.min(rowBuffer.length, MAX_BATCH));

  const { methodStr, minRep, maxRep, minIter, maxIter, approxBytes } = summarizeBatch(batch);
  const id = ++flushSeq;

  dbLog(`⤴︎ [mess #${id}] transmitting ${batch.length} rows (~${approxBytes}B) · methods: ${methodStr} · reps:[${minRep ?? "—"}-${maxRep ?? "—"}] · iters:[${minIter ?? "—"}-${maxIter ?? "—"}]`);

  inflight = true;
  const t0 = (typeof performance !== "undefined" ? performance.now() : Date.now());
  try {
    const result = await postJSON(ENDPOINTS.rows(jobId), { rows: batch }, { preferBeacon: !!final, keepalive: !!final });
    const t1 = (typeof performance !== "undefined" ? performance.now() : Date.now());
    const ms = Math.round(t1 - t0);
    retryDelay = RETRY_BASE_MS;
    dbLog(`[sent #${id}] ok via ${result.via} in ${ms}ms · remaining buffer: ${rowBuffer.length}`);
  } catch (err) {
    const { permanent, reason } = classifyPermanentRowsError(err);
    if (permanent) {
      rowsEndpointDisabled = true;
      rowsEndpointErrorReason = reason || "permanent server error";
      postMessage({ type: "log", line:
`🛑 rows endpoint permanently failing (${rowsEndpointErrorReason}); dumping current + future rows as artifacts.
   Check API_BASE="${API_BASE}" and ensure the backend exposes ${ENDPOINTS.rows("<jobId>")} and the required table exists.` });
      rowBuffer.unshift(...batch);
      const dump = rowBuffer.splice(0, Math.min(rowBuffer.length, MAX_BATCH));
      emitRowsAsArtifact(dump, rowsEndpointErrorReason);
    } else {
      rowBuffer.unshift(...batch);
      postMessage({ type: "log", line: `❌ [sent #${id}] ${String(err)} — retrying in ${retryDelay}ms (buffer: ${rowBuffer.length} rows)` });
      setTimeout(() => void flushNow(jobId, { final }), retryDelay);
      retryDelay = Math.min(Math.floor(retryDelay * 1.8), RETRY_MAX_MS);
    }
  } finally { inflight = false; }

  if (rowBuffer.length) setTimeout(() => void flushNow(jobId, { final }), 0);
  return true;
}

async function drainAll(jobId) {
  postMessage({ type: "log", line: `⤴︎ [final drain] start · buffer=${rowBuffer.length}` });
  while (inflight || rowBuffer.length) {
    if (!inflight && rowBuffer.length) {
      await flushNow(jobId, { final: true });
    } else {
      await new Promise(r => setTimeout(r, 50));
    }
  }
  postMessage({ type: "log", line: `[transmission complete]` });
}

/* ---------------- [EM_AllInfo] handling ---------------- */
function normalizeAAObject(raw) {
  if (!raw || typeof raw !== "object") return null;

  // Accept your C++ JSON keys; keep flexible for future renames.
  const out = {};

  // Required-ish fields
  const rep  = Number(raw.rep);
  const iter = Number(raw.iter);
  const root = raw.root != null ? String(raw.root) : "";

  if (Number.isFinite(rep)) out.rep = rep;
  if (Number.isFinite(iter)) out.iter = iter;
  out.root = root;

  // LLs
  if (raw.aa_ll_init  != null) out.aa_ll_init  = Number(raw.aa_ll_init);
  if (raw.aa_ll_final != null) out.aa_ll_final = Number(raw.aa_ll_final);

  // Arrays / matrices — pass through as-is (must be valid JSON-serializable)
  if (raw.root_prob_init  != null) out.aa_root_prob_init  = raw.root_prob_init;
  if (raw.root_prob_final != null) out.aa_root_prob_final = raw.root_prob_final;

  if (raw.exchangeability_matrix_init  != null) out.aa_exchangeability_init  = raw.exchangeability_matrix_init;
  if (raw.exchangeability_matrix_final != null) out.aa_exchangeability_final = raw.exchangeability_matrix_final;

  // Per-iter trace (map/object or array is fine)
  if (raw.aa_ll_per_iter != null) out.aa_ll_per_iter = raw.aa_ll_per_iter;

  // Keep the original payload for debugging
  out.raw_json = raw;

  return out;
}

async function uploadAAInfo(jobId, aaObj) {
  if (!jobId || !aaObj) return { ok: false };

  // Fill rep if missing from stream context
  if (!Number.isFinite(aaObj.rep) && Number.isFinite(lastRepSeen)) {
    aaObj.rep = lastRepSeen;
  }

  // Minimal validation: need job_id, root, rep to satisfy unique key (job_id, root, rep)
  const rootOk = typeof aaObj.root === "string" && aaObj.root.length > 0;
  const repOk  = Number.isFinite(aaObj.rep);
  if (!rootOk || !repOk) {
    postMessage({ type: "log", line: `[AA] skipped upload (root/rep missing). root=${aaObj.root ?? "(null)"} rep=${aaObj.rep ?? "(null)"}` });
    return { ok: false };
  }

  const payload = {
    job_id: jobId,
    rep: aaObj.rep,
    iter: Number.isFinite(aaObj.iter) ? aaObj.iter : null,
    root: aaObj.root || "",

    aa_ll_init:  Number.isFinite(aaObj.aa_ll_init)  ? aaObj.aa_ll_init  : null,
    aa_ll_final: Number.isFinite(aaObj.aa_ll_final) ? aaObj.aa_ll_final : null,

    aa_root_prob_init:  aaObj.aa_root_prob_init  ?? null,
    aa_root_prob_final: aaObj.aa_root_prob_final ?? null,

    aa_exchangeability_init:  aaObj.aa_exchangeability_init  ?? null,
    aa_exchangeability_final: aaObj.aa_exchangeability_final ?? null,

    aa_ll_per_iter: aaObj.aa_ll_per_iter ?? null,

    raw_json: aaObj.raw_json ?? null,
  };

  try {
    const size = approxBytesOfJson(payload);
    upLog(`↗ upsert /api/aa-info (job=${jobId}, root=${payload.root}, rep=${payload.rep}, ~${size}B)`);
    const res = await postJSON(ENDPOINTS.aa(), payload, { keepalive: true });
    dbLog(`[AA] saved via ${res.via}`);
    upLog(`✔ upsert /api/aa-info ok (via ${res.via})`);
    return { ok: true };
  } catch (e) {
    postMessage({ type: "log", line: `❌ /api/aa-info failed: ${String(e)}` });
    return { ok: false, error: e };
  }
}


function normalizeBestObject(raw) {
  if (!raw || typeof raw !== "object") return null;
  const out = {};
  if (raw.parsimony) out.parsimony = raw.parsimony;
  if (raw.dirichlet) out.dirichlet = raw.dirichlet;
  if (raw.hss)       out.hss       = raw.hss;
  if (!out.parsimony && raw.Parsimony) out.parsimony = raw.Parsimony;
  if (!out.dirichlet && raw.Dirichlet) out.dirichlet = raw.Dirichlet;
  if (!out.hss && raw.hss)             out.hss       = raw.hss;
  if (!out.parsimony && !out.dirichlet && !out.hss) {
    const m = String(raw.method || "").toLowerCase();
    if (m.includes("pars")) out.parsimony = raw;
    else if (m.includes("dir")) out.dirichlet = raw;
    else if (m.includes("hss")) out.hss = raw;
  }
  return (out.parsimony || out.dirichlet || out.hss) ? out : null;
}

function emitBestAsArtifact(bundle, reason) {
  try {
    const name = `em_allinfo_fallback_${jobIdForThisRun}_${++artifactSeq}.json`;
    const bytes = new TextEncoder().encode(JSON.stringify({ reason, job_id: jobIdForThisRun, bundle }));
    postMessage({ type: "log", line: `📦 Best-params saved as artifact "${name}" (${bytes.length}B) due to: ${reason}` });
    postMessage({ type: "artifact", name, bytes }, [bytes.buffer]);
  } catch (e) { postMessage({ type: "log", line: `⚠️ failed to emit best-params artifact: ${String(e)}` }); }
}

async function uploadBestMethods(jobId, shaped, { preferBeacon = false } = {}) {
  if (bestEndpointMissing) return { okAny: false, okAll: false };
  const methods = ["parsimony", "dirichlet", "hss"].filter(k => shaped[k]);
  if (!methods.length) return { okAny: false, okAll: false };

  let okAny = false, okAll = true;

  for (const k of methods) {
    const payload = { job_id: jobId, [k]: shaped[k] };
    const size = approxBytesOfJson(payload);
    upLog(`↗ upsert /api/allinfo (job=${jobId}, method=${k}, ~${size}B)`);
    dbLog(`→ [EM_AllInfo] uploading ${k} (~${size}B)`);
    try {
      const res = await postJSON(ENDPOINTS.allinfo(), payload, {
        preferBeacon: preferBeacon && size <= BEACON_MAX_BYTES,
        keepalive: !!preferBeacon
      });
      okAny = true;
      dbLog(`[EM_AllInfo] ${k} ok via ${res.via}`);
      upLog(`✔ upsert /api/allinfo (${k}) ok (via ${res.via})`);
    } catch (e) {
      okAll = false;
      if (e && (e.permanent || e.status === 404 || e.status === 410)) {
        bestEndpointMissing = true;
        postMessage({ type: "log", line: `[EM_AllInfo] endpoint missing (HTTP ${e.status}); suppressing further attempts.` });
        break;
      }
      postMessage({ type: "log", line: `❌ upsert /api/allinfo (${k}) failed: ${String(e)}` });
    }
  }

  return { okAny, okAll };
}

function attachRepIfMissing(obj) {
  if (!obj) return;
  const r = Number(obj.rep);
  if (!Number.isFinite(r) && Number.isFinite(lastRepSeen)) obj.rep = lastRepSeen;
}

async function flushPendingBestNow(reason = "manual flush") {
  if (pendingRep === null) return { okAny: false, okAll: false };
  const bundle = {};
  if (pendingBest.parsimony) bundle.parsimony = pendingBest.parsimony;
  if (pendingBest.dirichlet) bundle.dirichlet = pendingBest.dirichlet;
  if (pendingBest.hss)       bundle.hss       = pendingBest.hss;
  if (!Object.keys(bundle).length) return { okAny: false, okAll: false };

  lastBestObj = bundle;

  const result = await uploadBestMethods(jobIdForThisRun, bundle, { preferBeacon: false });
  bestUploaded = result.okAll;

  if (!result.okAny && !bestEndpointMissing) {
    emitBestAsArtifact({ rep: pendingRep, ...bundle }, `upload failed during ${reason}`);
  }

  pendingBest.parsimony = pendingBest.dirichlet = pendingBest.hss = null;
  pendingRep = null;

  return result;
}

function scheduleBestUpload() {
  if (bestFlushTimer) return;
  bestFlushTimer = setTimeout(() => {
    bestFlushTimer = null;
    void flushPendingBestNow("debounced flush");
  }, BEST_UPLOAD_DEBOUNCE_MS);
}

function handleBestLine(jsonText) {
  try {
    const shaped = normalizeBestObject(JSON.parse(jsonText));
    if (!shaped) return;

    attachRepIfMissing(shaped.parsimony);
    attachRepIfMissing(shaped.dirichlet);
    attachRepIfMissing(shaped.hss);

    const repNow = Number(
      (shaped.parsimony && shaped.parsimony.rep) ??
      (shaped.dirichlet && shaped.dirichlet.rep) ??
      (shaped.hss && shaped.hss.rep)
    );
    if (!Number.isFinite(repNow)) {
      postMessage({ type: "log", line: "[EM_AllInfo] missing rep; skipping upload" });
      return;
    }

    if (pendingRep !== null && repNow !== pendingRep) {
      void flushPendingBestNow("rep changed");
    }
    pendingRep = repNow;

    if (shaped.parsimony) pendingBest.parsimony = shaped.parsimony;
    if (shaped.dirichlet) pendingBest.dirichlet = shaped.dirichlet;
    if (shaped.hss)       pendingBest.hss       = shaped.hss;

    scheduleBestUpload();
  } catch (e) {
    postMessage({ type: "log", line: `EM_AllInfo-parse-error: ${String(e)} :: ${jsonText.slice(0,200)}` });
  }
}

/* ---------------- line printer ---------------- */

function handlePrint(txt) {
  const line = String(txt).trim();
  if (!line) return;

  if (jobIdForThisRun && line.startsWith("[ROW")) {
    const brace = line.indexOf("{");
    if (brace !== -1) {
      try { const obj = JSON.parse(line.slice(brace)); queueRow(jobIdForThisRun, obj); return; }
      catch (e) { postMessage({ type: "log", line: `row-parse-error: ${String(e)} :: ${line.slice(0,200)}` }); return; }
    }
  }

  if (jobIdForThisRun && line.startsWith("[EM_AllInfo]")) {
    const brace = line.indexOf("{");
    if (brace !== -1) { handleBestLine(line.slice(brace)); return; }
  }

  // NEW: handle AA optimization bundles
  if (jobIdForThisRun && line.startsWith("[AA_All_info]")) {
    const brace = line.indexOf("{");
    if (brace !== -1) {
      try {
        const raw = JSON.parse(line.slice(brace));
        const aaObj = normalizeAAObject(raw);
        if (aaObj) { void uploadAAInfo(jobIdForThisRun, aaObj); }
      } catch (e) {
        postMessage({ type: "log", line: `AA-parse-error: ${String(e)} :: ${line.slice(0,200)}` });
      }
      return;
    }
  }

  if (jobIdForThisRun) {
    const parsed = maybeParseTreeLine(line, jobIdForThisRun);
    if (parsed && parsed.newick) {
      const p = postTree(jobIdForThisRun, parsed.newick, parsed.meta || {});
      pendingTreeUploads.add(p);
      p.finally(() => pendingTreeUploads.delete(p));
    }
  }
  postMessage ({ type: "log", line: line});
}

/* ---------------- worker entry + runtime toggle handler ---------------- */

self.onmessage = async (e) => {
  // --- Toggles (unchanged) ---------------------------------------------------
  if (e?.data && e.data.__cmd === "setDbLogging") {
    DB_LOGS = !!e.data.enabled;
    try { console.debug(`[DB] logging ${DB_LOGS ? "ENABLED" : "DISABLED"}`); } catch {}
    return;
  }
  if (e?.data && e.data.__cmd === "setUpsertUiLogging") {
    UPSERT_UI_LOGS = !!e.data.enabled;
    postMessage({ type: "log", line: `ℹ︎ Upsert UI logs ${UPSERT_UI_LOGS ? "ENABLED" : "DISABLED"}` });
    return;
  }
  if (e?.data && e.data.__cmd === "setTreeUploading") {
    UPLOAD_TREES = !!e.data.enabled;
    postMessage({ type: "log", line: `Tree upload ${UPLOAD_TREES ? "ENABLED" : "DISABLED"}` });
    return;
  }
  if (e?.data && e.data.__cmd === "setApiSecret") {
    const s = String(e.data.secret || "");    
    API_SECRET = s;    
    AUTH = s && (/^Bearer\s+/i.test(s) ? s : `Bearer ${s}`);
    postMessage({ type: "log", line: API_SECRET ? "API secret set" : "API secret cleared" });
    return;
  }

  
  if (e?.data && e.data.__cmd === "runWithJson") {
  const { json, files, dnaSeqBytes, aaSeqBytes, jobId } = e.data;

  startedAtMs = Date.now();
  jobIdForThisRun = jobId || `job-${Date.now()}`;
  rowsEndpointDisabled = false;
  rowsEndpointErrorReason = "";
  treesUploadedCount = 0;

  // reset per-run state
  lastRepSeen = null;
  lastBestObj = null;
  bestUploaded = false;
  bestEndpointMissing = false;
  pendingRep = null;
  pendingBest.parsimony = pendingBest.dirichlet = pendingBest.hss = null;
  if (bestFlushTimer) { clearTimeout(bestFlushTimer); bestFlushTimer = null; }
  seenTreeHashesByJob.set(jobIdForThisRun, new Set());
  pendingTreeUploads.clear();

  try {
    const v = (self.NEXT_PUBLIC_COMMIT_SHA || Date.now());
    const createEmtr = (await import(`/wasm/emtr.js?v=${v}`)).default;

    const Module = await createEmtr({
      locateFile: (p) => `/wasm/${p}?v=${v}`,
      print: handlePrint,
      printErr: (txt) => { try { handlePrint(txt); } catch {} postMessage({ type: "log", line: `ERR: ${txt}` }); },
    });

    // ---- Parse config from the JSON so we can create the job row
    let cfg = null;
    try {
      const payload = JSON.parse(String(json || "{}"));
      const dna = payload?.settings?.dna;
      if (dna && typeof dna === "object") {
        cfg = {
          thr: Number(dna.thr),
          reps: Number(dna.reps),
          maxIter: Number(dna.maxIter),
          D_pi: Array.isArray(dna.D_pi) ? dna.D_pi.map(Number) : null,
          D_M:  Array.isArray(dna.D_M)  ? dna.D_M.map(Number)  : null,
        };
      }
    } catch {}

    // ---- Ensure the job exists BEFORE any PATCH/rows
    if (cfg && Number.isFinite(cfg.thr) && Number.isFinite(cfg.reps) && Number.isFinite(cfg.maxIter)) {
      await upsertJobMetadata(jobIdForThisRun, cfg);
    } else {
      postMessage({ type: "log", line: "⚠️ could not parse job config from JSON; skipping /api/jobs upsert" });
    }

    const FS = Module.FS;
    try { FS.mkdir("/work"); } catch {}
    try { FS.mkdir("/work/out"); } catch {}

    // Stage uploaded bytes (optional) to the EXACT paths referenced in JSON
    if (dnaSeqBytes && files?.dna) {
      FS.writeFile(files.dna, new Uint8Array(dnaSeqBytes));
      postMessage({ type: "log", line: `DNA → ${files.dna} (${(dnaSeqBytes.byteLength/1024).toFixed(1)} KiB)` });
    }
    if (aaSeqBytes && files?.aa) {
      FS.writeFile(files.aa, new Uint8Array(aaSeqBytes));
      postMessage({ type: "log", line: `AA  → ${files.aa} (${(aaSeqBytes.byteLength/1024).toFixed(1)} KiB)` });
    }

    // Line buffering (optional)
    if (typeof Module.set_line_buffered === "function") {
      Module.set_line_buffered();
    } else {
      try {
        const setLineBufferedC = Module.cwrap("set_line_buffered_c", "void", []);
        if (setLineBufferedC) setLineBufferedC();
      } catch {}
    }

    // Dispatch to the single-JSON C++ entrypoint
    const run_with_json = Module.cwrap("run_with_json", "number", ["string"]);
    if (!run_with_json) throw new Error("run_with_json not exported");
    postMessage({ type: "log", line: "WASM dispatch: run_with_json(...)" });

    // Use an allowed status value (backend rejects 'running')
    try {
      await setJobStatus(jobIdForThisRun, { status: "started", started_at: new Date().toISOString() });
    } catch {}

    const rc = run_with_json(json);

    // Drain row stream (if any), ship artifacts, flush "best", wait for trees
    try { await drainAll(jobIdForThisRun); } catch {}
    try {
      const outDir = "/work/out";
      const names = FS.readdir(outDir).filter((n) => n !== "." && n !== "..");
      for (const name of names) {
        const path = `${outDir}/${name}`;
        const st = FS.stat(path);
        if ((st.mode & 0o170000) === 0o100000) {
          const bytes = FS.readFile(path, { encoding: "binary" });
          postMessage({ type: "artifact", name, bytes }, [bytes.buffer]);

          if (/\.(newick|nwk|tree)$/i.test(name)) {
            try {
              const text = FS.readFile(path, { encoding: "utf8" });
              if (text && typeof text === "string" && text.includes("(") && text.includes(")")) {
                const p = postTree(jobIdForThisRun, text.trim(), { label: name });
                pendingTreeUploads.add(p); p.finally(() => pendingTreeUploads.delete(p));
              }
            } catch {}
          }
        }
      }
    } catch {}
    try { await flushPendingBestNow("end-of-run"); } catch {}
    if (lastBestObj && !bestUploaded && !bestEndpointMissing) {
      dbLog("↻ [EM_AllInfo] retrying at end-of-run…");
      try {
        const { okAny, okAll } = await uploadBestMethods(jobIdForThisRun, lastBestObj, { preferBeacon: true });
        if (!okAny) emitBestAsArtifact(lastBestObj, "network error while posting EM_AllInfo");
        bestUploaded = okAll;
      } catch {}
    }
    try { await awaitTreeUploads(); } catch {}

    const elapsedMs = Date.now() - startedAtMs;
    const minutes = (elapsedMs / 60000).toFixed(2);
    try {
      await setJobStatus(jobIdForThisRun, {
        status: rowsEndpointDisabled ? "blocked" : "completed",
        finished_at: new Date().toISOString(),
        dur_minutes: minutes,
        blocked_reason: rowsEndpointDisabled ? rowsEndpointErrorReason : undefined,
      });
    } catch {}
    postMessage({ type: "log", line: `⏱ finished in ${minutes}m` });
    postMessage({ type: "done", rc });
  } catch (err) {
    try { await drainAll(jobIdForThisRun); } catch {}
    try { await flushPendingBestNow("run failure"); } catch {}
    if (lastBestObj && !bestUploaded && !bestEndpointMissing) {
      dbLog("↻ [EM_AllInfo] retrying after failure…");
      try {
        const { okAny } = await uploadBestMethods(jobIdForThisRun, lastBestObj, { preferBeacon: true });
        if (!okAny) emitBestAsArtifact(lastBestObj, "network error while posting EM_AllInfo (after failure)");
      } catch {}
    }
    try { await awaitTreeUploads(); } catch {}

    const elapsedMs = Date.now() - startedAtMs;
    const minutes = (elapsedMs / 60000).toFixed(2);
    postMessage({ type: "log", line: `⏱ aborted after ${minutes}m` });
    postMessage({ type: "log", line: `❌ ${String(err?.message ?? err)}` });
    try {
      await setJobStatus(jobIdForThisRun, {
        status: rowsEndpointDisabled ? "blocked" : "failed",
        finished_at: new Date().toISOString(),
        dur_minutes: minutes,
        blocked_reason: rowsEndpointDisabled ? rowsEndpointErrorReason : undefined,
      });
    } catch {}
    postMessage({ type: "done", rc: 1 });
  }
  return;
}

  const { params, dnaSeqBytes, aaSeqBytes, jobId, fileNames } = e.data;

  startedAtMs = Date.now();
  jobIdForThisRun = jobId || `job-${Date.now()}`;
  rowsEndpointDisabled = false;
  rowsEndpointErrorReason = "";
  treesUploadedCount = 0;
  epsilonParam = Number.isFinite(Number(params?.epsilon)) ? Number(params.epsilon) : null;

  // reset row + best state
  lastRepSeen = null;
  lastBestObj = null;
  bestUploaded = false;
  bestEndpointMissing = false;
  pendingRep = null;
  pendingBest.parsimony = pendingBest.dirichlet = pendingBest.hss = null;
  if (bestFlushTimer) { clearTimeout(bestFlushTimer); bestFlushTimer = null; }

  // reset tree dedupe & pending trackers for this job
  seenTreeHashesByJob.set(jobIdForThisRun, new Set());
  pendingTreeUploads.clear();

  const v = (self.NEXT_PUBLIC_COMMIT_SHA || Date.now());
  const createEmtr = (await import(`/wasm/emtr.js?v=${v}`)).default;

  const Module = await createEmtr({
    locateFile: (p) => `/wasm/${p}?v=${v}`,
    print: handlePrint,
    printErr: (txt) => { try { handlePrint(txt); } catch {} postMessage({ type: "log", line: `ERR: ${txt}` }); },
  });

  async function shipArtifactsAndUploadTrees() {
    try {
      const outDir = "/work/out";
      const FS = Module.FS;
      const names = FS.readdir(outDir).filter((n) => n !== "." && n !== "..");
      for (const name of names) {
        const path = `${outDir}/${name}`;
        const st = FS.stat(path);
        if ((st.mode & 0o170000) === 0o100000) {
          const bytes = FS.readFile(path, { encoding: "binary" });
          postMessage({ type: "artifact", name, bytes }, [bytes.buffer]);

          if (/\.(newick|nwk|tree)$/i.test(name)) {
            try {
              const text = FS.readFile(path, { encoding: "utf8" });
              if (text && typeof text === "string" && text.includes("(") && text.includes(")")) {
                const p = postTree(jobIdForThisRun, text.trim(), { label: name });
                pendingTreeUploads.add(p); p.finally(() => pendingTreeUploads.delete(p));
              }
            } catch {}
          }
        }
      }
    } catch {}
  }

  try {
    const FS = Module.FS;
    try { FS.mkdir("/work"); } catch {}
    try { FS.mkdir("/work/out"); } catch {}

    // Stage inputs using ORIGINAL names (sanitized); pass empty strings if absent
    const dnaBase = sanitize(fileNames?.dna) || (dnaSeqBytes ? "dna_input.fa" : null);
    const aaBase  = sanitize(fileNames?.aa)  || (aaSeqBytes  ? "aa_input.faa" : null);

    let dnaPath = "", aaPath = "";
    if (dnaSeqBytes) {
      dnaPath = `/work/${dnaBase}`;
      FS.writeFile(dnaPath, new Uint8Array(dnaSeqBytes));
      postMessage({ type: "log", line: `DNA → ${dnaPath} (${(dnaSeqBytes.byteLength/1024).toFixed(1)} KiB)` });
    }
    if (aaSeqBytes) {
      aaPath = `/work/${aaBase}`;
      FS.writeFile(aaPath, new Uint8Array(aaSeqBytes));
      postMessage({ type: "log", line: `AA  → ${aaPath} (${(aaSeqBytes.byteLength/1024).toFixed(1)} KiB)` });
    }

    // Set line buffering (C or C++ bridge)
    if (typeof Module.set_line_buffered === "function") {
      Module.set_line_buffered();
    } else {
      try {
        const setLineBufferedC = Module.cwrap("set_line_buffered_c", "void", []);
        if (setLineBufferedC) setLineBufferedC();
      } catch {}
    }

    // Params
    const thr = Number(params.thr);
    const reps = params.reps | 0;
    const maxIter = params.maxIter | 0;
    const D_pi = Array.isArray(params?.D_pi) && params.D_pi.length === 4 ? params.D_pi.map(Number) : [100,100,100,100];
    const D_M  = Array.isArray(params?.D_M)  && params.D_M.length  === 4 ? params.D_M.map(Number)  : [100,2,2,2];

    await upsertJobMetadata(jobIdForThisRun, { thr, reps, maxIter, D_pi, D_M });

    postMessage({ type: "log", line: `WASM dispatch: EM_main (reps=${reps}, maxIter=${maxIter}, thr=${thr})` });

    // Run EM and stream rows
    let rc = 0;

    // Prefer embind class if present
    if (typeof Module.EMManager === "function") {
      const mgr = new Module.EMManager(
        dnaPath || "",        // dna path or empty
        aaPath  || "",        // aa path or empty
        reps,
        maxIter,
        thr,
        D_pi[0], D_pi[1], D_pi[2], D_pi[3],
        D_M[0],  D_M[1],  D_M[2],  D_M[3]
      );
      mgr.EM_main();
    } else {
      // C entrypoint fallback
      const EM_main_entry_c = Module.cwrap(
        "EM_main_entry_c",
        "number",
        ["string","string","number","number","number",
         "number","number","number","number",
         "number","number","number","number"]
      );
      if (!EM_main_entry_c) throw new Error("EM_main_entry_c not exported");
      rc = EM_main_entry_c(
        dnaPath || "",
        aaPath  || "",
        reps, maxIter, thr,
        D_pi[0], D_pi[1], D_pi[2], D_pi[3],
        D_M[0],  D_M[1],  D_M[2],  D_M[3]
      );
    }

    // Finish sending row-stream
    await drainAll(jobIdForThisRun);

    // Ship artifacts and attempt Newick uploads from files
    await shipArtifactsAndUploadTrees();

    // Flush any remaining best-bundle for the last rep
    await flushPendingBestNow("end-of-run");

    // Retry last bundle using Beacon if needed
    if (lastBestObj && !bestUploaded && !bestEndpointMissing) {
      dbLog("↻ [EM_AllInfo] retrying at end-of-run…");
      const { okAny, okAll } = await uploadBestMethods(jobIdForThisRun, lastBestObj, { preferBeacon: true });
      if (!okAny) emitBestAsArtifact(lastBestObj, "network error while posting EM_AllInfo");
      bestUploaded = okAll;
    }

    // Wait for in-flight tree uploads
    await awaitTreeUploads();

    const elapsedMs = Date.now() - startedAtMs;
    const minutes = (elapsedMs / 60000).toFixed(2);

    await setJobStatus(jobIdForThisRun, {
      status: rowsEndpointDisabled ? "blocked" : "completed",
      finished_at: new Date().toISOString(),
      dur_minutes: minutes,
      blocked_reason: rowsEndpointDisabled ? rowsEndpointErrorReason : undefined,
    });

    postMessage({ type: "log", line: `⏱ finished in ${minutes}m` });
    postMessage({ type: "done", rc });
  } catch (err) {
    // Finish sending row-stream even on failure
    await drainAll(jobIdForThisRun);

    // Flush whatever best bundle we have for the last rep
    await flushPendingBestNow("run failure");

    // Retry best bundle with Beacon once more
    if (lastBestObj && !bestUploaded && !bestEndpointMissing) {
      dbLog("↻ [EM_AllInfo] retrying after failure…");
      const { okAny } = await uploadBestMethods(jobIdForThisRun, lastBestObj, { preferBeacon: true });
      if (!okAny) emitBestAsArtifact(lastBestObj, "network error while posting EM_AllInfo (after failure)");
    }

    // Wait for trees even on failure
    await awaitTreeUploads();

    const elapsedMs = Date.now() - startedAtMs;
    const minutes = (elapsedMs / 60000).toFixed(2);
    postMessage({ type: "log", line: `⏱ aborted after ${minutes}m` });
    postMessage({ type: "log", line: `❌ ${String(err)}` });

    await setJobStatus(jobIdForThisRun, {
      status: rowsEndpointDisabled ? "blocked" : "failed",
      finished_at: new Date().toISOString(),
      dur_minutes: minutes,
      blocked_reason: rowsEndpointDisabled ? rowsEndpointErrorReason : undefined,
    });

    postMessage({ type: "done", rc: 1 });
  }
};
