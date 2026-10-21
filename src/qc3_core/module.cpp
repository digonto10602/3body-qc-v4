// qc3._core : compiled backend of the cache-free single-energy QC (prompts/001).
// Step 1: essentially empty module (version information only).
// Built by `./build.sh core` (alias `pyext`) into src/python/qc3/.
#include <pybind11/pybind11.h>

#ifndef QC3_GIT_COMMIT
#define QC3_GIT_COMMIT "unknown"
#endif

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
  m.doc() = "qc3 compiled core (step 1: empty module)";
  m.def("version", []() { return std::string("0.1.0-step1"); },
        "Version string of the compiled core.");
  m.def("git_commit", []() { return std::string(QC3_GIT_COMMIT); },
        "Short git commit the extension was built from ('unknown' if not in a git checkout).");
}
