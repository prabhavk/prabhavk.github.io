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

static std::string g_dayhoff_path;

// Core runner: constructs EMManager, applies Dayhoff path (if any), then runs EM
static int run_em_main_pipeline(const std::string& dnaSeqPath, const std::string& aaSeqPath,
                                int num_repetitions,
                                int max_iter,
                                double conv_threshold,
                                double pi_a_1, double pi_a_2, double pi_a_3, double pi_a_4,
                                double M_a_1, double M_a_2, double M_a_3, double M_a_4) {
  set_line_buffered_impl();
  EMManager mgr(dnaSeqPath, aaSeqPath, num_repetitions, max_iter, conv_threshold,
                pi_a_1, pi_a_2, pi_a_3, pi_a_4,
                M_a_1,  M_a_2,  M_a_3,  M_a_4);
  
  mgr.setDayhoffMatrixPath("/assets/Dayhoff.dat");
  mgr.EM_main();
  return 0;
}

extern "C" {

EMSCRIPTEN_KEEPALIVE
void set_line_buffered_c() { set_line_buffered_impl(); }

// Global C API to set the Dayhoff matrix path from JS:
//   Module.cwrap("setDayhoffMatrixPath","void",["string"])("/data/Dayhoff.dat");
EMSCRIPTEN_KEEPALIVE
void setDayhoffMatrixPath(const char* path) {
  g_dayhoff_path = path ? std::string(path) : std::string();
}

EMSCRIPTEN_KEEPALIVE
int EM_main_entry_c(const char* dnaSeqPath,
                    const char* aaSeqPath,
                    int num_repetitions,
                    int max_iter,
                    double conv_threshold,
                    double pi_a_1, double pi_a_2, double pi_a_3, double pi_a_4,
                    double M_a_1, double M_a_2, double M_a_3, double M_a_4) {
  return run_em_main_pipeline(std::string(dnaSeqPath ? dnaSeqPath : ""),
                              std::string(aaSeqPath  ? aaSeqPath  : ""),
                              num_repetitions, max_iter, conv_threshold,
                              pi_a_1, pi_a_2, pi_a_3, pi_a_4,
                              M_a_1,  M_a_2,  M_a_3,  M_a_4);
}

} // extern "C"

// ---------- Embind (class API only; no ambiguous helpers) ----------
// NOTE: We bind the method under the JS name "setDayhoffMatrixPath" even though the C++
// method name is SetRateMatrixPath, so your worker can call:
//    mgr.setDayhoffMatrixPath("/data/Dayhoff.dat");
EMSCRIPTEN_BINDINGS(emtr_module) {
  emscripten::class_<EMManager>("EMManager")
    // constructor(dnaPath, aaPath, reps, maxIter, thr, pi1..pi4, M1..M4)
    .constructor<std::string, std::string, int, int, double,
                 double,double,double,double,
                 double,double,double,double>()
    .function("EM_main",               &EMManager::EM_main)
    .function("EMparsimony",           &EMManager::EMparsimony)
    .function("EMdirichlet",           &EMManager::EMdirichlet)
    .function("EMssh",                 &EMManager::EMssh)
    .function("setDayhoffMatrixPath",  &EMManager::setDayhoffMatrixPath);
}
