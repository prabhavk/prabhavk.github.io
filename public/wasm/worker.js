// public/wasm/worker.js

// ---------- Upload config ----------
const API_BASE = "";                // "" => same origin; or e.g. "https://emtr-web.vercel.app"
const AUTH = "";                    // optional: "Bearer <token)" if you add auth later

const MAX_BATCH = 200;               // flush threshold (row count per POST)
const FLUSH_INTERVAL_MS = 1000;     // idle flush (debounce)
const RETRY_BASE_MS = 1500;         // backoff start
const RETRY_MAX_MS = 15000;         // backoff cap

// ---------- Upload state ----------
let jobIdForThisRun = "";
let rowBuffer = [];
let flushTimer = null;
let inflight = false;
let retryDelay = RETRY_BASE_MS;
let startedAtMs = 0;                // ⏱ timing
let flushSeq = 0;                   // sequential id for each flush

// ---------- helpers ----------
async function upsertJobMetadata(jobId, cfg) {
  try {
    const res = await fetch(`${API_BASE}/api/jobs`, {
      method: "POST",
      headers: {
        "content-type": "application/json",
        ...(AUTH ? { authorization: AUTH } : {}),
      },
      keepalive: true,
      body: JSON.stringify({
        job_id: jobId,
        thr: cfg.thr,
        reps: cfg.reps,
        max_iter: cfg.maxIter,
        ssh_rounds: cfg.sshRounds ?? cfg.ssh_rounds ?? 1, // ← NEW
        D_pi: cfg.D_pi,   // [a1,a2,a3,a4]
        D_M: cfg.D_M,     // [b1,b2,b3,b4]
        status: "started",
      }),
    });
    if (!res.ok) {
      const text = await res.text().catch(() => "");
      throw new Error(`HTTP ${res.status} ${res.statusText} :: ${text.slice(0, 300)}`);
    }
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
      postMessage({
        type: "log",
        line: `⚠️ job patch failed: HTTP ${res.status} :: ${text.slice(0,200)}`
      });
    } else {
      if (payload?.status)
        postMessage({ type: "log", line: `ℹ︎ job status -> ${payload.status}` });
    }
  } catch (e) {
    postMessage({ type: "log", line: `⚠️ job patch error: ${String(e)}` });
  }
}

// Approximate size of one row's JSON (for logging only)
function approxBytesOfRow(r) {
  return 8 + JSON.stringify(r).length + 1; // tag + json + newline
}

function summarizeBatch(batch) {
  const methods = Object.create(null);
  let minRep = Infinity, maxRep = -Infinity;
  let minIter = Infinity, maxIter = -Infinity;
  let approxBytes = 0;

  for (const r of batch) {
    const m = (r.method ?? "(none)").toString();
    methods[m] = (methods[m] || 0) + 1;

    if (Number.isFinite(r.rep)) {
      if (r.rep < minRep) minRep = r.rep;
      if (r.rep > maxRep) maxRep = r.rep;
    }
    if (Number.isFinite(r.iter)) {
      if (r.iter < minIter) minIter = r.iter;
      if (r.iter > maxIter) maxIter = r.iter;
    }
    approxBytes += approxBytesOfRow(r);
  }

  if (!Number.isFinite(minRep)) { minRep = null; maxRep = null; }
  if (!Number.isFinite(minIter)) { minIter = null; maxIter = null; }

  const methodStr = Object.keys(methods).length
    ? Object.entries(methods).map(([k,v]) => `${k}:${v}`).join(", ")
    : "—";

  return { methodStr, minRep, maxRep, minIter, maxIter, approxBytes };
}

// Coerce one row to canonical shape
function normalizeRow(raw) {
  const toNum = (v) => (v == null || v === "" ? NaN : Number(v));
  const toInt = (v) => (v == null || v === "" ? NaN : parseInt(v, 10));

  const ll_init_src =
    raw.ll_init ??
    raw.ll_initial ??
    raw.logLikelihood_initial ??
    raw.loglikelihood_initial;

  const ecd_first_src =
    raw.ecd_ll_first ??
    raw.logLikelihood_ecd_first ??
    raw.loglikelihood_ecd_first;

  const ecd_final_src =
    raw.ecd_ll_final ??
    raw.logLikelihood_ecd_final ??
    raw.loglikelihood_ecd_final;

  const ll_final_src =
    raw.ll_final ??
    raw.logLikelihood_final ??
    raw.loglikelihood_final;

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

// Flush a chunk (<= MAX_BATCH). Returns true if attempted a POST.
async function flushNow(jobId, { final }) {
  // if a flush is in-flight, let pending timer handle the next one
  if (inflight) return false;

  if (flushTimer) { clearTimeout(flushTimer); flushTimer = null; }
  if (!rowBuffer.length) return false;

  const batch = rowBuffer.splice(0, Math.min(rowBuffer.length, MAX_BATCH));
  const { methodStr, minRep, maxRep, minIter, maxIter, approxBytes } = summarizeBatch(batch);
  const id = ++flushSeq;

  postMessage({
    type: "log",
    line: `⤴︎ [#${id}] transmitting ${batch.length} rows (~${approxBytes}B) · methods: ${methodStr} · reps:[${minRep ?? "—"}-${maxRep ?? "—"}] · iters:[${minIter ?? "—"}-${maxIter ?? "—"}]`
  });

  inflight = true;
  const t0 = (typeof performance !== "undefined" ? performance.now() : Date.now());
  try {
    const res = await fetch(`${API_BASE}/api/jobs/${encodeURIComponent(jobId)}/rows`, {
      method: "POST",
      headers: {
        "content-type": "application/json",
        ...(AUTH ? { authorization: AUTH } : {}),
      },
      // keepalive can drop large bodies; use it only for the tiny final chunks
      keepalive: !!final,
      body: JSON.stringify({ rows: batch }),
    });

    const t1 = (typeof performance !== "undefined" ? performance.now() : Date.now());
    const ms = Math.round(t1 - t0);

    if (!res.ok) {
      const text = await res.text().catch(() => "");
      throw new Error(`HTTP ${res.status} ${res.statusText} :: ${text.slice(0, 400)}`);
    }

    retryDelay = RETRY_BASE_MS;
    postMessage({
      type: "log",
      line: `[#${id}] ok ${res.status} in ${ms}ms · remaining buffer: ${rowBuffer.length}`
    });
  } catch (err) {
    // put batch back and retry later
    rowBuffer.unshift(...batch);
    postMessage({
      type: "log",
      line: `❌ [flush #${id}] ${String(err)} — retrying in ${retryDelay}ms (buffer: ${rowBuffer.length} rows)`
    });
    setTimeout(() => void flushNow(jobId, { final }), retryDelay);
    retryDelay = Math.min(Math.floor(retryDelay * 1.8), RETRY_MAX_MS);
  } finally {
    inflight = false;
  }

  // If more rows are waiting, schedule another immediate flush
  if (rowBuffer.length) setTimeout(() => void flushNow(jobId, { final }), 0);
  return true;
}

// Drain everything at the end (keeps flushing until buffer empty & no inflight)
async function drainAll(jobId) {
  postMessage({ type: "log", line: `⤴︎ [final drain] start · buffer=${rowBuffer.length}` });
  while (inflight || rowBuffer.length) {
    if (!inflight && rowBuffer.length) {
      await flushNow(jobId, { final: true });
    } else {
      await new Promise(r => setTimeout(r, 50));
    }
  }
  postMessage({ type: "log", line: `[final drain] done` });
}

// Accept any [ROW...] tag; parse JSON from first '{'
function handlePrint(txt) {
  const line = String(txt).trim();
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
  postMessage({ type: "log", line });
}

// ---------- WASM worker entry ----------
self.onmessage = async (e) => {
  const { params, seqBytes, topoBytes, jobId } = e.data;
  startedAtMs = Date.now();                           // ⏱ start timing
  jobIdForThisRun = jobId || `job-${Date.now()}`;

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
    const sshRounds = (params.sshRounds | 0) || 1; // ← NEW
    const D_pi = Array.isArray(params?.D_pi) && params.D_pi.length === 4 ? params.D_pi.map(Number) : [100,100,100,100];
    const D_M  = Array.isArray(params?.D_M)  && params.D_M.length  === 4 ? params.D_M.map(Number)  : [100,2,2,2];

    // ① Save job config (status=started)
    await upsertJobMetadata(jobIdForThisRun, { thr, reps, maxIter, sshRounds, D_pi, D_M });

    postMessage({ type: "log", line: `WASM dispatch: EM_main (reps=${reps}, maxIter=${maxIter}, sshRounds=${sshRounds}, thr=${thr})` });

    // ② Run EM and stream rows
    let rc = 0;
    if (typeof Module.EMManager === "function") {
      // NEW ORDER: sequence, topology, reps, maxIter, ssh_rounds, thr, pi..., M...
      const mgr = new Module.EMManager(
        "/work/sequence.phyx",
        "/work/topology.csv",
        reps,
        maxIter,
        sshRounds,
        thr,
        D_pi[0], D_pi[1], D_pi[2], D_pi[3],
        D_M[0],  D_M[1],  D_M[2],  D_M[3]
      );
      mgr.EM_main();
    } else {
      // Keep the arity/types list; pass values in the NEW order.
      const EM_main_entry_c = Module.cwrap(
        "EM_main_entry_c",
        "number",
        [
          "string","string", // seq, topo
          "number","number","number","number", // reps, maxIter, ssh_rounds, thr
          "number","number","number","number", // pi1..pi4
          "number","number","number","number"  // M1..M4
        ]
      );
      if (!EM_main_entry_c) throw new Error("EM_main_entry_c not exported");
      rc = EM_main_entry_c(
        "/work/sequence.phyx",
        "/work/topology.csv",
        reps,
        maxIter,
        sshRounds,
        thr,
        D_pi[0], D_pi[1], D_pi[2], D_pi[3],
        D_M[0],  D_M[1],  D_M[2],  D_M[3]
      );
    }

    // Drain any remaining rows before finalization
    await drainAll(jobIdForThisRun);
    shipArtifacts();

    // ⏱ success timing
    const elapsedMs = Date.now() - startedAtMs;
    const minutes = (elapsedMs / 60000).toFixed(2);
    postMessage({ type: "log", line: `⏱ finished in ${minutes}m` });
    postMessage({ type: "done", rc });
    await setJobStatus(jobIdForThisRun, {
      status: "completed",
      finished_at: new Date().toISOString(),
      dur_minutes: minutes,
    });
  } catch (err) {
    // Ensure we try to drain on failure, too
    await drainAll(jobIdForThisRun);

    const elapsedMs = Date.now() - startedAtMs;
    const minutes = (elapsedMs / 60000).toFixed(2);
    postMessage({ type: "log", line: `⏱ aborted after ${minutes}m` });
    postMessage({ type: "log", line: `❌ ${String(err)}` });
    postMessage({ type: "done", rc: 1 });
    await setJobStatus(jobIdForThisRun, {
      status: "failed",
      finished_at: new Date().toISOString(),
      dur_minutes: minutes,
    });
  }
};
