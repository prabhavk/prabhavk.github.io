// public/wasm/worker.js

// (A) --- NEW: batching + API config ---------------------------------
const API_BASE = "";                // "" => same origin; or e.g. "http://localhost:3000" / "https://api.yourdomain.com"
const AUTH = "";                    // optional: bearer token if your API requires it

const rowBuffer = [];
let flushTimer = null;
let jobIdForThisRun;                // set from the page

function queueRow(jobId, row) {
  rowBuffer.push(row);
  if (rowBuffer.length >= 200) flushNow(jobId);
  if (!flushTimer) flushTimer = setTimeout(() => flushNow(jobId), 1000);
}

async function flushNow(jobId) {
  if (flushTimer) { clearTimeout(flushTimer); flushTimer = null; }
  const batch = rowBuffer.splice(0, rowBuffer.length);
  if (!batch.length) return;
  try {
    await fetch(`${API_BASE}/api/jobs/${jobId}/rows`, {
      method: "POST",
      headers: {
        "content-type": "application/json",
        ...(AUTH ? { authorization: `Bearer ${AUTH}` } : {}),
      },
      body: JSON.stringify({ rows: batch }),
    });
    postMessage({ type: "log", line: `⤴︎ uploaded ${batch.length} rows` });
  } catch (err) {
    postMessage({ type: "log", line: `upload failed: ${String(err)} — retrying soon` });
    rowBuffer.unshift(...batch);
    if (!flushTimer) flushTimer = setTimeout(() => flushNow(jobId), 2000);
  }
}

// --------------------------------------------------------------------

self.onmessage = async (e) => {
  const { params, seqBytes, topoBytes, jobId } = e.data;   // (B) NEW: accept jobId
  jobIdForThisRun = jobId;

  // Cache-bust in dev so you don't run stale wasm
  const v = Date.now();
  const createEmtr = (await import(`/wasm/emtr.js?v=${v}`)).default;

  const Module = await createEmtr({
    locateFile: (p) => `/wasm/${p}?v=${v}`,
    // (C) NEW: parse emitted rows in stdout and buffer them
    print: (txt) => {
      const line = String(txt).trim();
      if (jobIdForThisRun && line.startsWith("[ROW]\t")) {
        try {
          const row = JSON.parse(line.slice(6)); // after "[ROW]\t"
          queueRow(jobIdForThisRun, row);
          return; // don't spam UI log with data rows
        } catch (e) {
          postMessage({ type: "log", line: `row-parse-error: ${String(e)} :: ${line}` });
        }
      }
      postMessage({ type: "log", line: txt });
    },
    printErr: (txt) => postMessage({ type: "log", line: `ERR: ${txt}` }),
  });

  // Helper to ship any artifacts from /work/out
  function shipArtifacts() {
    try {
      const outDir = "/work/out";
      const FS = Module.FS;
      const names = FS.readdir(outDir).filter((n) => n !== "." && n !== "..");
      for (const name of names) {
        const path = `${outDir}/${name}`;
        const stat = FS.stat(path);
        // Only send regular files
        if ((stat.mode & 0o170000) === 0o100000) {
          const bytes = FS.readFile(path, { encoding: "binary" }); // Uint8Array
          postMessage({ type: "artifact", name, bytes }, [bytes.buffer]);
        }
      }
    } catch {
      /* no-op if /work/out missing or empty */
    }
  }

  try {
    const FS = Module.FS;
    try { FS.mkdir("/work"); } catch {}
    try { FS.mkdir("/work/out"); } catch {}

    // Write inputs where your bridge & embind code expect them
    FS.writeFile("/work/sequence.phyx", new Uint8Array(seqBytes));
    FS.writeFile("/work/topology.csv", new Uint8Array(topoBytes));

    // Set line-buffered stdout/stderr for progressive logs
    if (typeof Module.set_line_buffered === "function") {
      Module.set_line_buffered(); // embind version
    } else {
      const setLineBufferedC = Module.cwrap("set_line_buffered_c", "void", []);
      if (setLineBufferedC) setLineBufferedC();
    }

    const thr = Number(params.thr);
    const reps = params.reps | 0;
    const maxIter = params.maxIter | 0;
    const fmt = params.seqFormat || "phylip";

    postMessage({ type: "log", line: `WASM dispatch: method="${params.method}"` });

    let rc = 0;

    if (params.method === "main") {
      // Prefer the embind class path (no-arg method) if available
      if (typeof Module.EMManager === "function") {
        const mgr = new Module.EMManager(
          "/work/sequence.phyx", fmt, "/work/topology.csv", "/work/out",
          reps, maxIter, thr
        );
        mgr.EM_main(); // no args
        await flushNow(jobIdForThisRun);  // (D) NEW: final flush before finishing
        shipArtifacts();
        postMessage({ type: "done", rc: 0 });
        return;
      }

      // Fallback to the C-export bridge
      const EM_main_entry_c = Module.cwrap(
        "EM_main_entry_c",
        "number",
        ["string", "string", "number", "number", "number", "string"]
      );
      if (!EM_main_entry_c) {
        throw new Error("EM_main_entry_c not exported");
      }
      rc = EM_main_entry_c("/work/sequence.phyx", "/work/topology.csv", thr, reps, maxIter, fmt);
      await flushNow(jobIdForThisRun);    // (D) NEW: final flush
      shipArtifacts();
      postMessage({ type: "done", rc });
      return;
    }

    // Non-main methods: embind class API
    if (typeof Module.EMManager === "function") {
      const mgr = new Module.EMManager(
        "/work/sequence.phyx", fmt, "/work/topology.csv", "/work/out",
        reps, maxIter, thr
      );
      switch (params.method) {
        case "parsimony": mgr.EMparsimony(); break;
        case "ssh":       mgr.EMssh();       break;
        case "dirichlet":
        default:          mgr.EMdirichlet(); break;
      }
      await flushNow(jobIdForThisRun);      // (D) NEW: final flush
      shipArtifacts();
      postMessage({ type: "done", rc: 0 });
      return;
    }

    // If we got here, nothing was callable
    throw new Error("No EMManager binding; and method !== 'main' for bridge fallback.");
  } catch (err) {
    // try to flush whatever we buffered before signaling error
    await flushNow(jobIdForThisRun);        // (D) NEW: flush on error too
    postMessage({ type: "log", line: `❌ ${String(err)}` });
    postMessage({ type: "done", rc: 1 });
  }
};
