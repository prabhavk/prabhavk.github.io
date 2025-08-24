// public/wasm/worker.js

// ---------- Upload config ----------
const API_BASE = self.API_BASE || ""; // "" => same origin; or e.g. "https://emtr-web.vercel.app"
const AUTH = "";                      // optional: "Bearer <token)"

// Allow overriding endpoint paths without rebuilding the worker
const ENDPOINTS = {
  rows: (jobId) => `${API_BASE}/api/jobs/${encodeURIComponent(jobId)}/rows`,
  allinfo: () => `${API_BASE}${self.ALL_INFO_PATH || "/api/allinfo"}`,
};

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

// Track [EM_AllInfo] handling
let bestSeenCount = 0;            // how many times we successfully sent
let lastBestObj = null;           // last parsed best object
let bestUploaded = false;         // whether last best object was uploaded OK

// Coalesce + suppress for [EM_AllInfo]
let pendingBestObj = null;
let bestFlushTimer = null;
const BEST_UPLOAD_DEBOUNCE_MS = 500;
let bestEndpointMissing = false;  // set true after a 404/410 so we stop trying

// If the rows endpoint is permanently failing (e.g., table missing), stop retrying
let rowsEndpointDisabled = false;
let rowsEndpointErrorReason = "";
let artifactSeq = 0;

// ---------- helpers ----------
function approxBytesOfRow(r) {
  // crude but stable estimate
  return 8 + JSON.stringify(r).length + 1;
}

function approxBytesOfJson(obj) {
  return JSON.stringify(obj).length;
}

function canUseBeacon() {
  try {
    return typeof navigator !== "undefined" && typeof navigator.sendBeacon === "function";
  } catch {
    return false;
  }
}

function tryBeacon(url, obj) {
  try {
    const blob = new Blob([JSON.stringify(obj)], { type: "application/json" });
    return navigator.sendBeacon(url, blob);
  } catch {
    return false;
  }
}

async function postJSON(url, obj, { preferBeacon = false, keepalive = false } = {}) {
  // Small payload + allowed? try Beacon first
  if (preferBeacon && canUseBeacon() && approxBytesOfJson(obj) <= BEACON_MAX_BYTES) {
    const ok = tryBeacon(url, obj);
    if (ok) return { ok: true, via: "beacon" };
    // fall through to fetch
  }

  // Fallback to fetch (optionally keepalive)
  const res = await fetch(url, {
    method: "POST",
    headers: {
      "content-type": "application/json",
      ...(AUTH ? { authorization: AUTH } : {}),
    },
    body: JSON.stringify(obj),
    keepalive,
    cache: "no-store",
  });
  if (!res.ok) {
    const text = await res.text().catch(() => "");
    const err = new Error(`HTTP ${res.status} ${res.statusText} :: ${text.slice(0, 400)}`);
    err.status = res.status;
    err.bodyText = text || "";
    err.permanent = res.status === 404 || res.status === 410;
    throw err;
  }
  return { ok: true, via: "fetch" };
}

function classifyPermanentRowsError(err) {
  if (!err) return { permanent: false, reason: "" };
  const status = err.status | 0;
  const msg = (err.bodyText || err.message || "").toLowerCase();

  // 404/410 are permanent for sure
  if (status === 404 || status === 410) {
    return { permanent: true, reason: `endpoint returned ${status}` };
  }

  // Common "table missing" signals -> permanent until server is fixed
  const tableMissing =
    msg.includes("doesn't exist") ||
    msg.includes("does not exist") ||
    msg.includes("no such table") ||
    msg.includes("unknown table") ||
    msg.includes("table") && msg.includes("not found");

  if (status === 500 && tableMissing) {
    return { permanent: true, reason: "target table is missing on the server" };
  }

  return { permanent: false, reason: "" };
}

async function upsertJobMetadata(jobId, cfg) {
  try {
    await postJSON(`${API_BASE}/api/jobs`, {
      job_id: jobId,
      thr: cfg.thr,
      reps: cfg.reps,
      max_iter: cfg.maxIter,
      D_pi: cfg.D_pi,
      D_M: cfg.D_M,
      status: "started",
    }, { keepalive: true }); // harmless outside teardown
    postMessage({ type: "log", line: "ℹ︎ job config saved" });
  } catch (err) {
    postMessage({ type: "log", line: `⚠️ failed to save job config: ${String(err)}` });
  }
}

async function setJobStatus(jobId, payload) {
  try {
    const res = await fetch(`${API_BASE}/api/jobs/${encodeURIComponent(jobId)}`, {
      method: "PATCH",
      headers: {
        "content-type": "application/json",
        ...(AUTH ? { authorization: AUTH } : {}),
      },
      keepalive: true,
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const text = await res.text().catch(() => "");
      postMessage({ type: "log", line: `⚠️ job patch failed: HTTP ${res.status} :: ${text.slice(0,200)}` });
    } else {
      if (payload?.status) postMessage({ type: "log", line: `ℹ︎ job status -> ${payload.status}` });
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

  const ll_init_src =
    raw.ll_init ?? raw.ll_initial ?? raw.logLikelihood_initial ?? raw.loglikelihood_initial;
  const ecd_first_src =
    raw.ecd_ll_first ?? raw.logLikelihood_ecd_first ?? raw.loglikelihood_ecd_first;
  const ecd_final_src =
    raw.ecd_ll_final ?? raw.logLikelihood_ecd_final ?? raw.loglikelihood_ecd_final;
  const ll_final_src =
    raw.ll_final ?? raw.logLikelihood_final ?? raw.loglikelihood_final;

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
  if (!flushTimer) {
    flushTimer = setTimeout(() => void flushNow(jobIdForThisRun, { final: false }), FLUSH_INTERVAL_MS);
  }
}

function queueRow(jobId, row) {
  const r = normalizeRow(row);
  if (!r) return;
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
    if (batch.length > 0 && bytes + cost > maxBytes) break; // keep at least one row
    batch.push(rowBuffer.shift());
    bytes += cost;
  }
  return batch;
}

function encodeCsv(rows) {
  const header = [
    "job_id","method","root","rep","iter",
    "ll_init","ecd_ll_first","ecd_ll_final","ll_final"
  ];
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
    postMessage({
      type: "log",
      line: `📦 Rows saved as artifact "${name}" (${bytes.length}B) due to: ${reason}`
    });
    postMessage({ type: "artifact", name, bytes }, [bytes.buffer]);
  } catch (e) {
    postMessage({ type: "log", line: `⚠️ failed to emit fallback artifact: ${String(e)}` });
  }
}

async function flushNow(jobId, { final }) {
  if (inflight) return false;
  if (flushTimer) { clearTimeout(flushTimer); flushTimer = null; }
  if (!rowBuffer.length) return false;

  // If rows endpoint is disabled, just dump to artifacts and skip network
  if (rowsEndpointDisabled) {
    const batch = rowBuffer.splice(0, Math.min(rowBuffer.length, MAX_BATCH));
    const id = ++flushSeq;
    emitRowsAsArtifact(batch, rowsEndpointErrorReason || "rows endpoint disabled");
    postMessage({ type: "log", line: `[sent #${id}] skipped network (endpoint disabled) · remaining buffer: ${rowBuffer.length}` });
    if (rowBuffer.length) setTimeout(() => void flushNow(jobId, { final }), 0);
    return true;
  }

  // For final flush, keep the payload small enough for Beacon.
  const batch = final && canUseBeacon()
    ? takeByteCappedBatch(MAX_BATCH, BEACON_MAX_BYTES - 512) // small header margin
    : rowBuffer.splice(0, Math.min(rowBuffer.length, MAX_BATCH));

  const { methodStr, minRep, maxRep, minIter, maxIter, approxBytes } = summarizeBatch(batch);
  const id = ++flushSeq;

  postMessage({
    type: "log",
    line: `⤴︎ [mess #${id}] transmitting ${batch.length} rows (~${approxBytes}B) · methods: ${methodStr} · reps:[${minRep ?? "—"}-${maxRep ?? "—"}] · iters:[${minIter ?? "—"}-${maxIter ?? "—"}]`
  });

  inflight = true;
  const t0 = (typeof performance !== "undefined" ? performance.now() : Date.now());
  try {
    const result = await postJSON(
      ENDPOINTS.rows(jobId),
      { rows: batch },
      {
        // prefer Beacon only for "final" drains
        preferBeacon: !!final,
        keepalive: !!final, // helps if Beacon not available
      }
    );

    const t1 = (typeof performance !== "undefined" ? performance.now() : Date.now());
    const ms = Math.round(t1 - t0);
    retryDelay = RETRY_BASE_MS;
    postMessage({
      type: "log",
      line: `[sent #${id}] ok via ${result.via} in ${ms}ms · remaining buffer: ${rowBuffer.length}`
    });
  } catch (err) {
    const { permanent, reason } = classifyPermanentRowsError(err);
    if (permanent) {
      // Stop retrying network; route future batches to artifacts
      rowsEndpointDisabled = true;
      rowsEndpointErrorReason = reason || "permanent server error";
      postMessage({
        type: "log",
        line:
`🛑 rows endpoint permanently failing (${rowsEndpointErrorReason}); dumping current + future rows as artifacts.
   Check API_BASE="${API_BASE}" and ensure the backend exposes ${ENDPOINTS.rows("<jobId>")} and the required table exists.`,
      });
      // Put rows back then dump them immediately as artifact to preserve order
      rowBuffer.unshift(...batch);
      const dump = rowBuffer.splice(0, Math.min(rowBuffer.length, MAX_BATCH));
      emitRowsAsArtifact(dump, rowsEndpointErrorReason);
    } else {
      // transient: put rows back at the front (preserve ordering) and retry
      rowBuffer.unshift(...batch);
      postMessage({
        type: "log",
        line: `❌ [sent #${id}] ${String(err)} — retrying in ${retryDelay}ms (buffer: ${rowBuffer.length} rows)`
      });
      setTimeout(() => void flushNow(jobId, { final }), retryDelay);
      retryDelay = Math.min(Math.floor(retryDelay * 1.8), RETRY_MAX_MS);
    }
  } finally {
    inflight = false;
  }

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

function normalizeBestObject(raw) {
  // Expect shape like: { parsimony: {...}, dirichlet: {...}, ssh: {...} }
  if (!raw || typeof raw !== "object") return null;
  const out = {};

  // Pass-through preferred keys
  if (raw.parsimony) out.parsimony = raw.parsimony;
  if (raw.dirichlet) out.dirichlet = raw.dirichlet;
  if (raw.ssh)       out.ssh       = raw.ssh;

  // Accept a few alternates just in case
  if (!out.parsimony && raw.Parsimony) out.parsimony = raw.Parsimony;
  if (!out.dirichlet && raw.Dirichlet) out.dirichlet = raw.Dirichlet;
  if (!out.ssh && raw.SSH)             out.ssh       = raw.SSH;

  // Single-struct fallback
  if (!out.parsimony && !out.dirichlet && !out.ssh) {
    const method = String(raw.method || "").toLowerCase();
    if (method.includes("pars")) out.parsimony = raw;
    else if (method.includes("dir")) out.dirichlet = raw;
    else if (method.includes("ssh")) out.ssh = raw;
  }

  if (!out.parsimony && !out.dirichlet && !out.ssh) return null;
  return out;
}

async function uploadBest(jobId, bestRaw, { preferBeacon = false } = {}) {
  if (bestEndpointMissing) return false;

  const shaped = normalizeBestObject(bestRaw);
  if (!shaped) {
    postMessage({ type: "log", line: "⚠️ [EM_AllInfo] payload shape not recognized; skipping" });
    return false;
  }

  const payload = { job_id: jobId, ...shaped };
  const size = approxBytesOfJson(payload);
  postMessage({ type: "log", line: `→ [EM_AllInfo] upload start (~${size}B)` });

  try {
    const res = await postJSON(ENDPOINTS.allinfo(), payload, {
      preferBeacon,
      keepalive: preferBeacon, // helpful fallback if Beacon rejected
    });
    bestSeenCount += 1;
    postMessage({ type: "log", line: `✅ [EM_AllInfo] upserted (${bestSeenCount}) via ${res.via}` });
    return true;
  } catch (e) {
    if (e && (e.permanent || e.status === 404 || e.status === 410)) {
      bestEndpointMissing = true;
      postMessage({
        type: "log",
        line:
`🛑 [EM_AllInfo] endpoint not found (HTTP ${e.status}). Further attempts suppressed.
   Check API_BASE="${API_BASE}" and ALL_INFO_PATH="${self.ALL_INFO_PATH || "/api/allinfo"}" point to a server exposing that route.`,
      });
      return false;
    }
    postMessage({ type: "log", line: `❌ [EM_AllInfo] upload error: ${String(e)}` });
    return false;
  }
}

function scheduleBestUpload() {
  if (bestFlushTimer) return;
  bestFlushTimer = setTimeout(() => {
    const obj = pendingBestObj;
    pendingBestObj = null;
    bestFlushTimer = null;
    if (obj) {
      uploadBest(jobIdForThisRun, obj, { preferBeacon: false }).then((ok) => {
        bestUploaded = ok;
        if (!ok && !bestEndpointMissing) postMessage({ type: "log", line: "⏳ [EM_AllInfo] will retry at end" });
      });
    }
  }, BEST_UPLOAD_DEBOUNCE_MS);
}

function handlePrint(txt) {
  const line = String(txt).trim();
  if (!line) return;

  // Streamed rows
  if (jobIdForThisRun && line.startsWith("[ROW")) {
    const brace = line.indexOf("{");
    if (brace !== -1) {
      try {
        const parsed = JSON.parse(line.slice(brace));
        queueRow(jobIdForThisRun, parsed);
        return;
      } catch (e) {
        postMessage({ type: "log", line: `row-parse-error: ${String(e)} :: ${line.slice(0,200)}` });
      }
    }
  }

  // Best-parameters JSON
  if (jobIdForThisRun && line.startsWith("[EM_AllInfo]")) {
    const brace = line.indexOf("{");
    if (brace !== -1) {
      try {
        const parsed = JSON.parse(line.slice(brace));
        lastBestObj = parsed;
        bestUploaded = false;

        // Coalesce mid-run updates to avoid spamming the server
        pendingBestObj = parsed;
        scheduleBestUpload();

        postMessage({ type: "log", line: "… [EM_AllInfo] parsed" });
        return;
      } catch (e) {
        postMessage({ type: "log", line: `EM_AllInfo-parse-error: ${String(e)} :: ${line.slice(0,200)}` });
      }
    }
  }

  // Everything else to UI log
  postMessage({ type: "log", line });
}

// ---------- WASM worker entry ----------
self.onmessage = async (e) => {
  const { params, seqBytes, topoBytes, jobId } = e.data;
  startedAtMs = Date.now();
  jobIdForThisRun = jobId || `job-${Date.now()}`;
  bestSeenCount = 0;
  lastBestObj = null;
  bestUploaded = false;
  bestEndpointMissing = false;
  rowsEndpointDisabled = false;
  rowsEndpointErrorReason = "";
  pendingBestObj = null;
  if (bestFlushTimer) { clearTimeout(bestFlushTimer); bestFlushTimer = null; }

  const v = (self.NEXT_PUBLIC_COMMIT_SHA || Date.now());
  const createEmtr = (await import(`/wasm/emtr.js?v=${v}`)).default;

  const Module = await createEmtr({
    locateFile: (p) => `/wasm/${p}?v=${v}`,
    print: handlePrint,
    printErr: (txt) => postMessage({ type: "log", line: `ERR: ${txt}` }),
  });

  function shipArtifacts() {
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
        }
      }
    } catch {}
  }

  try {
    const FS = Module.FS;
    try { FS.mkdir("/work"); } catch {}
    try { FS.mkdir("/work/out"); } catch {}

    FS.writeFile("/work/sequence.phyx", new Uint8Array(seqBytes));
    FS.writeFile("/work/topology.csv", new Uint8Array(topoBytes));

    if (typeof Module.set_line_buffered === "function") {
      Module.set_line_buffered();
    } else {
      const setLineBufferedC = Module.cwrap("set_line_buffered_c", "void", []);
      if (setLineBufferedC) setLineBufferedC();
    }

    const thr = Number(params.thr);
    const reps = params.reps | 0;
    const maxIter = params.maxIter | 0;
    const D_pi = Array.isArray(params?.D_pi) && params.D_pi.length === 4 ? params.D_pi.map(Number) : [100,100,100,100];
    const D_M  = Array.isArray(params?.D_M)  && params.D_M.length  === 4 ? params.D_M.map(Number)  : [100,2,2,2];

    await upsertJobMetadata(jobIdForThisRun, { thr, reps, maxIter, D_pi, D_M });

    postMessage({ type: "log", line: `WASM dispatch: EM_main (reps=${reps}, maxIter=${maxIter}, thr=${thr})` });

    // Run EM and stream rows
    let rc = 0;
    if (typeof Module.EMManager === "function") {
      const mgr = new Module.EMManager(
        "/work/sequence.phyx",
        "/work/topology.csv",
        reps,
        maxIter,
        thr,
        D_pi[0], D_pi[1], D_pi[2], D_pi[3],
        D_M[0],  D_M[1],  D_M[2],  D_M[3]
      );
      mgr.EM_main();
    } else {
      const EM_main_entry_c = Module.cwrap(
        "EM_main_entry_c",
        "number",
        ["string","string","number","number","number",
         "number","number","number","number",
         "number","number","number","number"]
      );
      if (!EM_main_entry_c) throw new Error("EM_main_entry_c not exported");
      rc = EM_main_entry_c(
        "/work/sequence.phyx", "/work/topology.csv",
        thr, reps, maxIter,
        D_pi[0], D_pi[1], D_pi[2], D_pi[3],
        D_M[0],  D_M[1],  D_M[2],  D_M[3]
      );
    }

    await drainAll(jobIdForThisRun);
    shipArtifacts();

    // Retry [EM_AllInfo] once at end — prefer Beacon here
    if (lastBestObj && !bestUploaded && !bestEndpointMissing) {
      postMessage({ type: "log", line: "↻ [EM_AllInfo] retrying at end-of-run…" });
      bestUploaded = await uploadBest(jobIdForThisRun, lastBestObj, { preferBeacon: true });
    }

    const elapsedMs = Date.now() - startedAtMs;
    const minutes = (elapsedMs / 60000).toFixed(2);
    postMessage({ type: "log", line: `⏱ finished in ${minutes}m` });
    postMessage({ type: "done", rc });
    await setJobStatus(jobIdForThisRun, {
      status: rowsEndpointDisabled ? "blocked" : "completed",
      finished_at: new Date().toISOString(),
      dur_minutes: minutes,
      blocked_reason: rowsEndpointDisabled ? rowsEndpointErrorReason : undefined,
    });
  } catch (err) {
    await drainAll(jobIdForThisRun);

    // Retry [EM_AllInfo] once on failure, too (Beacon preferred)
    if (lastBestObj && !bestUploaded && !bestEndpointMissing) {
      postMessage({ type: "log", line: "↻ [EM_AllInfo] retrying after failure…" });
      bestUploaded = await uploadBest(jobIdForThisRun, lastBestObj, { preferBeacon: true });
    }

    const elapsedMs = Date.now() - startedAtMs;
    const minutes = (elapsedMs / 60000).toFixed(2);
    postMessage({ type: "log", line: `⏱ aborted after ${minutes}m` });
    postMessage({ type: "log", line: `❌ ${String(err)}` });
    postMessage({ type: "done", rc: 1 });
    await setJobStatus(jobIdForThisRun, {
      status: rowsEndpointDisabled ? "blocked" : "failed",
      finished_at: new Date().toISOString(),
      dur_minutes: minutes,
      blocked_reason: rowsEndpointDisabled ? rowsEndpointErrorReason : undefined,
    });
  }
};
