self.onmessage = async (e) => {
  const { params, seqBytes, topoBytes } = e.data;

  // load the Emscripten module the build produced
  const createEmtr = (await import('/wasm/emtr.js')).default;
  const mod = await createEmtr({
    locateFile: (p) => `/wasm/${p}`,
    print: (txt) => self.postMessage({ type: 'log', line: txt }),
    printErr: (txt) => self.postMessage({ type: 'log', line: `ERR: ${txt}` }),
  });

  try {
    // write inputs to the virtual FS
    mod.FS.writeFile('/sequence.phyx', new Uint8Array(seqBytes));
    mod.FS.writeFile('/topology.csv', new Uint8Array(topoBytes));

    // construct and run
    const mgr = new mod.EMManager(
      '/sequence.phyx',
      params.seqFormat || 'phylip',
      '/topology.csv',
      '/out',
      params.reps | 0,
      params.maxIter | 0,
      Number(params.thr)
    );

    if (params.method === 'parsimony')      mgr.EMparsimony();
    else if (params.method === 'ssh')       mgr.EMssh();
    else                                    mgr.EMdirichlet();


    self.postMessage({ type: 'done', rc: 0 });
  } catch (err) {
    self.postMessage({ type: 'log', line: `❌ ${String(err)}` });
    self.postMessage({ type: 'done', rc: 1 });
  }
};
