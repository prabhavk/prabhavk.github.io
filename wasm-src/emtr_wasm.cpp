// emtr_wasm.cpp (drop-in replacement)
#include "emtr_core.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include "third_party/json/single_include/nlohmann/json.hpp"

using json = nlohmann::json;
using std::string;
using std::vector;

// ---------- shared helpers ----------
static void set_line_buffered_impl() {
  setvbuf(stdout, nullptr, _IOLBF, 0);
  setvbuf(stderr, nullptr, _IOLBF, 0);
}

static std::string g_dayhoff_path;

static inline bool has_path(const json& j, const char* a) {
  return j.contains(a) && !j[a].is_null();
}
static inline bool has_path(const json& j, const char* a, const char* b) {
  return j.contains(a) && !j[a].is_null() && j[a].contains(b) && !j[a][b].is_null();
}
static inline string json_string_or(const json& j, const char* key, const string& fallback) {
  if (!has_path(j, key)) return fallback;
  return j.at(key).get<string>();
}
static inline string json_nested_string_or(const json& j, const char* a, const char* b, const string& fallback) {
  if (!has_path(j, a, b)) return fallback;
  return j.at(a).at(b).get<string>();
}

// Pull "settings" object if present; otherwise fallback to the root
static inline const json& root_settings_or_self(const json& cfg) {
  if (cfg.contains("settings") && cfg["settings"].is_object()) return cfg["settings"];
  return cfg;
}

// File path resolution compatible with both new and old payloads:
// - new: cfg.files.{dna,aa,dayhoff}
// - old: cfg.fileNames.{dna,aa} and optionally cfg.dayhoffPath
static inline void resolve_file_paths(const json& cfg, string& dnaPath, string& aaPath, string& dayhoffPath) {
  // defaults
  dnaPath     = "/data/RAxML_DNA_test.fa";
  aaPath      = "/data/RAxML_AA_test.fa";
  dayhoffPath = "/data/Dayhoff.dat";

  if (has_path(cfg, "files")) {
    if (has_path(cfg, "files", "dna"))     dnaPath     = cfg["files"]["dna"].get<string>();
    if (has_path(cfg, "files", "aa"))      aaPath      = cfg["files"]["aa"].get<string>();
    if (has_path(cfg, "files", "dayhoff")) dayhoffPath = cfg["files"]["dayhoff"].get<string>();
  }

  // Back-compat fallbacks
  if (dnaPath == "/data/RAxML_DNA_test.fa")
    dnaPath = json_nested_string_or(cfg, "fileNames", "dna", dnaPath);
  if (aaPath == "/data/RAxML_AA_test.fa")
    aaPath  = json_nested_string_or(cfg, "fileNames", "aa",  aaPath);
  if (dayhoffPath == "/data/Dayhoff.dat" && has_path(cfg, "dayhoffPath"))
    dayhoffPath = cfg["dayhoffPath"].get<string>();
}

// Build a minimal JSON from legacy C entry params to reuse the same runner
static inline string build_json_for_legacy_c(
    const char* dnaPath, const char* aaPath,
    int reps, int maxIter, double thr,
    double d_pi1, double d_pi2, double d_pi3, double d_pi4,
    double d_m1,  double d_m2,  double d_m3,  double d_m4) {

  json j;
  j["job_id"] = "legacy-c";
  j["files"] = {
    {"dna", dnaPath ? dnaPath : ""},
    {"aa",  aaPath  ? aaPath  : ""}
  };
  json S; // settings
  // DNA (EMDNA) block
  S["dna"] = {
    {"reps", reps},
    {"maxIter", maxIter},
    {"thrPerRegime", { {"coarse", thr}, {"medium", thr}, {"fine", thr} }},
    {"D_pi", { d_pi1, d_pi2, d_pi3, d_pi4 }},
    {"D_M",  { d_m1,  d_m2,  d_m3,  d_m4  }},
    {"requireAlphaGe1", true},
    {"includeParsimony", true},
    {"includeHSS", true}
  };
  // AA block left default unless needed by your core
  S["aa"] = json::object();
  j["settings"] = S;
  return j.dump();
}

// Core runner: constructs EMManager from parsed JSON
static int run_em_from_json(const std::string& settings_json) {
  set_line_buffered_impl();

  json cfg = json::parse(settings_json); // may throw
  const json& S = root_settings_or_self(cfg);

  string dnaPath, aaPath, dayhoffPath;
  resolve_file_paths(cfg, dnaPath, aaPath, dayhoffPath);

  // Example: pull a few common settings if ever needed
  int reps    = 10;
  int maxIter = 100;
  double thr  = 0.01;

  if (has_path(S, "dna")) {
    const auto& D = S["dna"];
    if (D.contains("reps"))     reps    = D["reps"].get<int>();
    if (D.contains("maxIter"))  maxIter = D["maxIter"].get<int>();
    if (has_path(D, "thrPerRegime", "medium"))
      thr = D["thrPerRegime"]["medium"].get<double>();
  }

  // Use your constructor that accepts (dnaPath, aaPath, settings_json)
  EMManager mgr(dnaPath, aaPath, settings_json);

  // Dayhoff path override: global > cfg.files.dayhoff > cfg.dayhoffPath > default
  if (!g_dayhoff_path.empty()) {
    mgr.setDayhoffMatrixPath(g_dayhoff_path);
  } else {
    mgr.setDayhoffMatrixPath(dayhoffPath);
  }

  // Kick off
  mgr.EM_main();
  return 0;
}

extern "C" {

// --- tiny utility ABI --------------------------------------------------------
EMSCRIPTEN_KEEPALIVE
void set_line_buffered_c() { set_line_buffered_impl(); }

EMSCRIPTEN_KEEPALIVE
void setDayhoffMatrixPath(const char* path) {
  g_dayhoff_path = path ? string(path) : string();
}

// --- NEW: canonical single-JSON entrypoint -----------------------------------
EMSCRIPTEN_KEEPALIVE
int run_with_json(const char* json_str) {
  if (!json_str) return -1;
  try {
    return run_em_from_json(string(json_str));
  } catch (const std::exception& e) {
    std::fprintf(stderr, "run_with_json failed: %s\n", e.what());
    return -2;
  } catch (...) {
    std::fprintf(stderr, "Unknown error in run_with_json\n");
    return -3;
  }
}

// --- Back-compat: legacy C entry takes explicit args -------------------------
EMSCRIPTEN_KEEPALIVE
int EM_main_entry_c(const char* dnaPath,
                    const char* aaPath,
                    int reps, int maxIter, double thr,
                    double d_pi1, double d_pi2, double d_pi3, double d_pi4,
                    double d_m1,  double d_m2,  double d_m3,  double d_m4) {
  try {
    const string j = build_json_for_legacy_c(
      dnaPath, aaPath, reps, maxIter, thr,
      d_pi1, d_pi2, d_pi3, d_pi4,
      d_m1,  d_m2,  d_m3,  d_m4
    );
    return run_em_from_json(j);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EM_main_entry_c failed: %s\n", e.what());
    return -2;
  } catch (...) {
    std::fprintf(stderr, "Unknown error in EM_main_entry_c\n");
    return -3;
  }
}

} // extern "C"

// ---------- Embind (class API only; no ambiguous helpers) ----------
EMSCRIPTEN_BINDINGS(emtr_module) {
  emscripten::class_<EMManager>("EMManager")
    .constructor<std::string, std::string, std::string>()
    .function("EM_main",               &EMManager::EM_main)
    .function("EMparsimony",           &EMManager::EMparsimony)
    .function("EMdirichlet",           &EMManager::EMdirichlet)
    .function("EMhss",                 &EMManager::EMhss);    
}
