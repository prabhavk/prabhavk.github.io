// public/wasm/worker.js

// ---------- Upload config ----------
const API_BASE = "";                // "" => same origin; or e.g. "https://emtr-web.vercel.app"
const AUTH = "";                    // optional: "Bearer <token)" if you add auth later

const MAX_BATCH = 200;              // flush threshold
const FLUSH_INTERVAL_MS = 1000;     // idle flush
const RETRY_BASE_MS = 1500;         // backoff start
const RETRY_MAX_MS = 15_000;        // backoff cap

// ---------- Upload state ----------
let jobIdForThisRun = "";
let rowBuffer = [];
let flushTimer = null;
let inflight = false;
let retryDelay = RETRY_BASE_MS;

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

async function setJobStatus(jobId, status) {
  try {
    await fetch(`${API_BASE}/api/jobs/${encodeURIComponent(jobId)}/status`, {
      method: "POST",
      headers: {
        "content-type": "application/json",
        ...(AUTH ? { authorization: AUTH } : {}),
      },
      keepalive: true, // helps if tab closes
      body: JSON.stringify({ status }),
    });
  } catch (e) {
    postMessage({ type: "log", line: `job status update failed: ${String(e)}` });
  }
}

// Coerce one row to canonical shape
function normalizeRow(raw) {
  const toNum = (v) => (v == null || v === "" ? NaN : Number(v));
  const toInt = (v) => (v == null || v === "" ? NaN : parseInt(v, 10));

  const method_src = raw.method;

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
    method: method_src,
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
    flushTimer = setTimeout(() => flushNow(jobIdForThisRun), FLUSH_INTERVAL_MS);
  }
}

function queueRow(jobId, row) {
  const r = normalizeRow(row);
  if (!r) return;
  rowBuffer.push(r);

  if (rowBuffer.length >= MAX_BATCH) {
    if (!inflight) flushNow(jobId);
    else scheduleFlush();
  } else {
    scheduleFlush();
  }
}

async function flushNow(jobId) {
  if (flushTimer) { clearTimeout(flushTimer); flushTimer = null; }
  if (inflight) return;
  if (!rowBuffer.length) return;

  const batch = rowBuffer.splice(0, rowBuffer.length);
  inflight = true;
  try {
    const res = await fetch(`${API_BASE}/api/jobs/${encodeURIComponent(jobId)}/rows`, {
      method: "POST",
      headers: {
        "content-type": "application/json",
        ...(AUTH ? { authorization: AUTH } : {}),
      },
      keepalive: true,
      body: JSON.stringify({ rows: batch }),
    });

    if (!res.ok) {
      const text = await res.text().catch(() => "");
      throw new Error(`HTTP ${res.status} ${res.statusText} :: ${text.slice(0, 400)}`);
    }

    retryDelay = RETRY_BASE_MS;
    postMessage({ type: "log", line: `⤴︎ uploaded ${batch.length} rows (remaining: ${rowBuffer.length})` });
  } catch (err) {
    rowBuffer.unshift(...batch);
    postMessage({ type: "log", line: `upload failed: ${String(err)} — retrying in ${retryDelay}ms` });
    setTimeout(() => flushNow(jobId), retryDelay);
    retryDelay = Math.min(Math.floor(retryDelay * 1.8), RETRY_MAX_MS);
  } finally {
    inflight = false;
  }

  // keep draining if more piled up
  if (rowBuffer.length) setTimeout(() => flushNow(jobId), 0);
}

// ---------- WASM worker entry ----------
self.onmessage = async (e) => {
  const { params, seqBytes, topoBytes, jobId } = e.data;
  jobIdForThisRun = jobId || `job-${Date.now()}`;

  const v = (self.NEXT_PUBLIC_COMMIT_SHA || Date.now());
  const createEmtr = (await import(`/wasm/emtr.js?v=${v}`)).default;

  const Module = await createEmtr({
    locateFile: (p) => `/wasm/${p}?v=${v}`,
    print: (txt) => {
      const line = String(txt).trim();
      if (jobIdForThisRun && line.startsWith("[ROW]\t")) {
        try {
          const parsed = JSON.parse(line.slice(6));
          queueRow(jobIdForThisRun, parsed);
          return;
        } catch (e) {
          postMessage({ type: "log", line: `row-parse-error: ${String(e)} :: ${line}` });
        }
      }
      postMessage({ type: "log", line: txt });
    },
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

    // ① Save job config (status=started)
    await upsertJobMetadata(jobIdForThisRun, { thr, reps, maxIter, D_pi, D_M });

    postMessage({ type: "log", line: `WASM dispatch: EM_main` });

    // ② Run EM and stream rows
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

    await flushNow(jobIdForThisRun);
    shipArtifacts();
    postMessage({ type: "done", rc });

    // ③ mark completed
    setJobStatus(jobIdForThisRun, "completed");
  } catch (err) {
    await flushNow(jobIdForThisRun);
    postMessage({ type: "log", line: `❌ ${String(err)}` });
    postMessage({ type: "done", rc: 1 });

    // ③ mark failed
    setJobStatus(jobIdForThisRun, "failed");
  }
};
