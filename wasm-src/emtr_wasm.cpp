#include "emtr_core.hpp"
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>

extern "C" {
EMSCRIPTEN_KEEPALIVE
void set_line_buffered() {
  setvbuf(stdout, nullptr, _IOLBF, 0);
  setvbuf(stderr, nullptr, _IOLBF, 0);
}
}

using namespace emscripten;

EMSCRIPTEN_BINDINGS(emtr_module) {
  class_<EMManager>("EMManager")
    .constructor<std::string, std::string, std::string, std::string, int, int, double>()
    .function("EM_main", &EMManager::EM_main)
    .function("EMparsimony", &EMManager::EMparsimony)
    .function("EMdirichlet", &EMManager::EMdirichlet)
    .function("EMssh", &EMManager::EMssh);
}
