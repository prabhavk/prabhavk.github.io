// emtr_wasm.cpp
#include "emtr_core.hpp"
#include <cstdio>
#include <string>
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>

// ---------- shared helpers ----------
static void set_line_buffered_impl() {
  setvbuf(stdout, nullptr, _IOLBF, 0);
  setvbuf(stderr, nullptr, _IOLBF, 0);
}

static int run_em_main_pipeline(const std::string& seqPath,
                                const std::string& topoPath,
                                double thr,
                                int reps,
                                int maxIter,
                                const std::string& seqFormat,
                                const std::string& prefix = "/work/out") {
  set_line_buffered_impl();
  // ctor: (sequence_file, seq_file_format, topology_file, prefix, num_reps, max_iter, conv_threshold)
  EMManager mgr(seqPath, seqFormat, topoPath, prefix, reps, maxIter, thr);
  mgr.EM_main(); // no-arg method using the state above
  return 0;
}

// ---------- C exports (usable via cwrap) ----------
extern "C" {

EMSCRIPTEN_KEEPALIVE
void set_line_buffered_c() { set_line_buffered_impl(); }

EMSCRIPTEN_KEEPALIVE
int EM_main_entry_c(const char* seqPath,
                    const char* topoPath,
                    double thr,
                    int reps,
                    int maxIter,
                    const char* seqFormat) {
  return run_em_main_pipeline(seqPath, topoPath, thr, reps, maxIter, seqFormat, "/work/out");
}

} // extern "C"

// ---------- Embind (class API + convenience fns) ----------
EMSCRIPTEN_BINDINGS(emtr_module) {
  emscripten::class_<EMManager>("EMManager")
    // (sequence_file, seq_file_format, topology_file, prefix, num_reps, max_iter, conv_threshold)
    .constructor<std::string, std::string, std::string, std::string, int, int, double>()
    .function("EM_main",     &EMManager::EM_main)
    .function("EMparsimony", &EMManager::EMparsimony)
    .function("EMdirichlet", &EMManager::EMdirichlet)
    .function("EMssh",       &EMManager::EMssh);

  // Embind-callable helpers:
  emscripten::function("set_line_buffered", &set_line_buffered_impl);

  emscripten::function(
    "EM_main_entry",
    emscripten::optional_override([](const std::string& seqPath,
                                     const std::string& topoPath,
                                     double thr,
                                     int reps,
                                     int maxIter,
                                     const std::string& seqFormat) {
      return run_em_main_pipeline(seqPath, topoPath, thr, reps, maxIter, seqFormat, "/work/out");
    })
  );
}
