// public/wasm/worker.js

// ---------- Upload config ----------
const API_BASE = "";                // "" => same origin; or e.g. "https://emtr-web.vercel.app"
const AUTH = "";                    // optional: "Bearer <token>" if you add auth later

const MAX_BATCH = 200;              // flush when buffer hits this size
const FLUSH_INTERVAL_MS = 1000;     // or after 1s of inactivity
const RETRY_BASE_MS = 1500;         // backoff start
const RETRY_MAX_MS = 15_000;        // backoff cap

// ---------- Upload state ----------
let jobIdForThisRun = "";
let rowBuffer = [];
let flushTimer = null;
let inflight = false;
let retryDelay = RETRY_BASE_MS;

// Coerce/validate a parsed row so it satisfies the API schema
function normalizeRow(raw) {
  const r = { ...raw };

  // If your API still requires "method", lock it to "main"
  r.method = "main";

  // ints required by the API
  r.rep = Number.isInteger(r.rep) ? r.rep : parseInt(r.rep ?? 0, 10);
  r.iter = Number.isInteger(r.iter) ? r.iter : parseInt(r.iter ?? 0, 10);

  // numeric fields
  r.ll_pars       = Number(r.ll_pars);
  r.edc_ll_first  = Number(r.edc_ll_first);
  r.edc_ll_final  = Number(r.edc_ll_final);
  r.ll_final      = Number(r.ll_final);

  return r;
}

function scheduleFlush() {
  if (!flushTimer) {
    flushTimer = setTimeout(() => flushNow(jobIdForThisRun), FLUSH_INTERVAL_MS);
  }
}

function queueRow(jobId, row) {
  rowBuffer.push(normalizeRow(row));
  if (rowBuffer.length >= MAX_BATCH) {
    flushNow(jobId);
  } else {
    scheduleFlush();
  }
}

async function flushNow(jobId) {
  if (flushTimer) { clearTimeout(flushTimer); flushTimer = null; }
  if (inflight) return;                 // don't overlap
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

    // success: reset backoff and notify UI
    retryDelay = RETRY_BASE_MS;
    postMessage({ type: "log", line: `⤴︎ uploaded ${batch.length} rows` });
  } catch (err) {
    // put rows back at the front of the buffer
    rowBuffer.unshift(...batch);
    postMessage({ type: "log", line: `upload failed: ${String(err)} — retrying in ${retryDelay}ms` });
    setTimeout(() => flushNow(jobId), retryDelay);
    retryDelay = Math.min(Math.floor(retryDelay * 1.8), RETRY_MAX_MS);
  } finally {
    inflight = false;
  }
}

// ---------- WASM worker entry ----------
self.onmessage = async (e) => {
  const { params, seqBytes, topoBytes, jobId } = e.data;
  jobIdForThisRun = jobId || `job-${Date.now()}`;

  // Cache-bust so we never load stale wasm in dev/preview
  const v = (self.NEXT_PUBLIC_COMMIT_SHA || Date.now());
  const createEmtr = (await import(`/wasm/emtr.js?v=${v}`)).default;

  const Module = await createEmtr({
    locateFile: (p) => `/wasm/${p}?v=${v}`,
    print: (txt) => {
      const line = String(txt).trim();
      // Expect lines like: [ROW]\t{...json...}
      if (jobIdForThisRun && line.startsWith("[ROW]\t")) {
        try {
          const parsed = JSON.parse(line.slice(6));
          queueRow(jobIdForThisRun, parsed);
          return; // don't echo data rows to UI
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
    } catch { /* ignore */ }
  }

  try {
    const FS = Module.FS;
    try { FS.mkdir("/work"); } catch {}
    try { FS.mkdir("/work/out"); } catch {}

    // Write inputs where the native side expects them
    FS.writeFile("/work/sequence.phyx", new Uint8Array(seqBytes));
    FS.writeFile("/work/topology.csv", new Uint8Array(topoBytes));

    // Line-buffered stdout for progressive logs
    if (typeof Module.set_line_buffered === "function") {
      Module.set_line_buffered();
    } else {
      const setLineBufferedC = Module.cwrap("set_line_buffered_c", "void", []);
      if (setLineBufferedC) setLineBufferedC();
    }

    const thr = Number(params.thr);
    const reps = params.reps | 0;
    const maxIter = params.maxIter | 0;

    // Fixed format (sequence format field removed from UI)
    const fmt = "phylip";

    postMessage({ type: "log", line: `WASM dispatch: EM_main` });

    // --- Single unified path: EM_main only ---
    let rc = 0;
    if (typeof Module.EMManager === "function") {
      const mgr = new Module.EMManager("/work/sequence.phyx", fmt, "/work/topology.csv", "/work/out", reps, maxIter, thr);
      mgr.EM_main(); // internally: parsimony → dirichlet → ssh
    } else {
      const EM_main_entry_c = Module.cwrap(
        "EM_main_entry_c",
        "number",
        ["string", "string", "number", "number", "number", "string"]
      );
      if (!EM_main_entry_c) throw new Error("EM_main_entry_c not exported");
      rc = EM_main_entry_c("/work/sequence.phyx", "/work/topology.csv", thr, reps, maxIter, fmt);
    }

    // Finalize: flush rows, ship artifacts, and finish
    await flushNow(jobIdForThisRun);
    shipArtifacts();
    postMessage({ type: "done", rc });
  } catch (err) {
    // Flush whatever we have before failing
    await flushNow(jobIdForThisRun);
    postMessage({ type: "log", line: `❌ ${String(err)}` });
    postMessage({ type: "done", rc: 1 });
  }
};
