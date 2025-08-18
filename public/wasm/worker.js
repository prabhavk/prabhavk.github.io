// public/wasm/worker.js
self.onmessage = async (e) => {
  const { params, seqBytes, topoBytes } = e.data;

  // Cache-bust in dev so you don't run stale wasm
  const v = Date.now();
  const createEmtr = (await import(`/wasm/emtr.js?v=${v}`)).default;

  const Module = await createEmtr({
    locateFile: (p) => `/wasm/${p}?v=${v}`,
    print: (txt) => postMessage({ type: "log", line: txt }),
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
      shipArtifacts();
      postMessage({ type: "done", rc: 0 });
      return;
    }

    // If we got here, nothing was callable
    throw new Error("No EMManager binding; and method !== 'main' for bridge fallback.");
  } catch (err) {
    postMessage({ type: "log", line: `❌ ${String(err)}` });
    postMessage({ type: "done", rc: 1 });
  }
};
