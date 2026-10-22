// qc3._core : compiled backend of the cache-free single-energy QC (prompts/001).
// Step 1: version info.  Step 2: oracle_point (the existing CPU path, unchanged physics, one energy, one core).
// Built by `./build.sh core` (alias `pyext`) into src/python/qc3/.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>

#include "oracle_point.hpp"
#include "fast_point.hpp"

#ifndef QC3_GIT_COMMIT
#define QC3_GIT_COMMIT "unknown"
#endif

namespace py = pybind11;

namespace {

// Eigen (column-major) -> C-contiguous row-major complex128 ndarray (copy).
py::array_t<std::complex<double>> to_numpy(const Eigen::MatrixXcd& M) {
  py::array_t<std::complex<double>> a({py::ssize_t(M.rows()), py::ssize_t(M.cols())});
  auto r = a.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < M.rows(); ++i)
    for (py::ssize_t j = 0; j < M.cols(); ++j) r(i, j) = M(i, j);
  return a;
}

qc3_oracle::Input parse_input(const py::dict& p) {
  qc3_oracle::Input in;
  in.irrep = p["irrep"].cast<std::string>();
  in.Lbyas = p["Lbyas"].cast<double>();
  in.xi = p["xi"].cast<double>();
  in.atmK = p["atmK"].cast<double>();
  in.atmpi = p["atmpi"].cast<double>();
  in.eta_1 = p["eta_1"].cast<double>();
  in.eta_2 = p["eta_2"].cast<double>();
  in.alpha = p["alpha"].cast<double>();
  in.epsilon_H = p["epsilon_H"].cast<double>();
  in.max_shell_num = p["max_shell_num"].cast<int>();
  in.waves_1 = p["waves_1"].cast<std::vector<int>>();
  in.waves_2 = p["waves_2"].cast<std::vector<int>>();
  in.parity = p["parity"].cast<int>();
  in.eig_tol = p["eig_tol"].cast<double>();
  in.a0_Kpi = p["a0_Kpi"].cast<double>();
  in.a1_Kpi = p["a1_Kpi"].cast<double>();
  in.a0_KK = p["a0_KK"].cast<double>();
  in.k2_pwave_form = p["k2_pwave_form"].cast<std::string>();
  in.Ecm = p["Ecm"].cast<double>();
  const auto c = p["c"].cast<std::vector<double>>();
  if (c.size() != 4) throw std::invalid_argument("oracle_point: c must have 4 entries");
  for (int i = 0; i < 4; ++i) in.c[std::size_t(i)] = c[std::size_t(i)];
  in.compute_k3 = p["compute_k3"].cast<bool>();
  return in;
}

py::dict output_to_dict(const qc3_oracle::Output& o) {

  py::dict d;
  d["Ecm"] = o.Ecm; d["E_lab"] = o.E_lab;
  d["N1"] = o.N1; d["N2"] = o.N2; d["N"] = o.N; d["Nproj"] = o.Nproj;
  d["nnP"] = py::make_tuple(o.nnP[0], o.nnP[1], o.nnP[2]);
  py::array_t<int> basis({py::ssize_t(o.basis.size()), py::ssize_t(6)});
  {
    auto b = basis.mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < py::ssize_t(o.basis.size()); ++i)
      for (int j = 0; j < 6; ++j) b(i, j) = o.basis[std::size_t(i)][std::size_t(j)];
  }
  d["basis"] = basis;   // columns: flavour, nx, ny, nz, ell, m
  d["F"] = to_numpy(o.F); d["G"] = to_numpy(o.G); d["K2inv"] = to_numpy(o.K2inv);
  d["F3"] = to_numpy(o.F3); d["F3inv"] = to_numpy(o.F3inv); d["Vsel"] = to_numpy(o.Vsel);
  d["K3"] = to_numpy(o.K3); d["Q"] = to_numpy(o.Q);
  if (o.has_k3_basis) {
    py::dict kb;
    const char* names[4] = {"iso0", "iso1", "B", "E"};
    for (int i = 0; i < 4; ++i) kb[names[i]] = to_numpy(o.K3_basis[std::size_t(i)]);
    d["K3_basis"] = kb;
  } else {
    d["K3_basis"] = py::none();
  }
  d["sign"] = o.sign; d["sign_imag"] = o.sign_imag; d["lnabsdet"] = o.lnabsdet;
  d["timings"] = o.timings;
  d["stats"] = o.stats;
  return d;
}

py::dict oracle_point(const py::dict& p) {
  const qc3_oracle::Input in = parse_input(p);
  qc3_oracle::Output o;
  {
    py::gil_scoped_release release;
    o = qc3_oracle::evaluate(in);
  }
  return output_to_dict(o);
}

py::dict fast_point(const py::dict& p, double trunc_eps, const std::string& projector, const std::vector<std::string>& want, bool g_symmetric) {
  const qc3_oracle::Input in = parse_input(p);
  qc3_fast::Options opt;
  opt.trunc_eps = trunc_eps;
  opt.projector = projector;
  opt.g_symmetric = g_symmetric;
  const auto has = [&](const char* name) { return std::find(want.begin(), want.end(), std::string(name)) != want.end(); };
  opt.want_F = has("F"); opt.want_G = has("G"); opt.want_K2inv = has("K2inv"); opt.want_F3 = has("F3");
  opt.want_F3inv = has("F3inv"); opt.want_Vsel = has("Vsel"); opt.want_K3 = has("K3") || in.compute_k3;
  qc3_oracle::Output o;
  {
    py::gil_scoped_release release;
    o = qc3_fast::evaluate(in, opt);
  }
  return output_to_dict(o);
}

py::list dimension_jumps(const py::dict& p, double E_lo, double E_hi) {
  const qc3_oracle::Input in = parse_input(p);
  std::vector<double> v;
  {
    py::gil_scoped_release release;
    v = qc3_fast::dimension_jumps(in, E_lo, E_hi);
  }
  py::list out;
  for (double x : v) out.append(x);
  return out;
}

// Test hook: the oracle's Kin (config_maker_4_momentum_first) vs the threshold builder's Kin, compared bit for bit.
py::dict compare_kinematics(const py::dict& p) {
  const qc3_oracle::Input in = parse_input(p);
  qc3_oracle::Kin a, b;
  {
    py::gil_scoped_release release;
    a = qc3_oracle::make_kinematics(in);
    b = qc3_fast::fast_kinematics(in);
  }
  bool same = (a.N1 == b.N1 && a.N2 == b.N2 && a.N == b.N && a.En_c == b.En_c && a.L == b.L && a.basis == b.basis &&
               a.np_config == b.np_config && a.nk_config == b.nk_config && a.plm_config == b.plm_config &&
               a.klm_config == b.klm_config && a.total_P == b.total_P);
  py::dict d;
  d["identical"] = same;
  d["N"] = a.N;
  d["N_fast"] = b.N;
  return d;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "qc3 compiled core (oracle_point; fast path to follow)";
  m.def("version", []() { return std::string("0.5.0-step7.6"); },
        "Version string of the compiled core.");
  m.def("git_commit", []() { return std::string(QC3_GIT_COMMIT); },
        "Short git commit the extension was built from ('unknown' if not in a git checkout).");
  m.def("clear_caches", []() { qc3_fast::clear_fast_caches(); qc3_fast::clear_threshold_tables(); },
        "Clear the fast path's in-memory caches (projector V_sel per complete key; basis threshold tables).");
  m.def("cache_info", []() {
    auto& c = qc3_fast::projector_cache();
    std::lock_guard<std::mutex> lock(c.mtx);
    py::dict d;
    d["projector_entries"] = (long long)c.map.size();
    d["projector_hits"] = c.hits;
    d["projector_misses"] = c.misses;
    d["threshold_entries"] = qc3_fast::threshold_entries();
    return d;
  }, "Counters of the fast path's in-memory caches.");
  m.def("fast_point", &fast_point, py::arg("params"), py::arg("trunc_eps") = 1e-16, py::arg("projector") = std::string("orbit"),
        py::arg("want") = std::vector<std::string>{"F", "G", "K2inv", "F3", "F3inv", "Vsel", "K3", "Q"}, py::arg("g_symmetric") = false,
        "Fast path (prompts/001 step 6, prompts/002 7.4-7.6); same params as oracle_point plus the truncation tolerance, the "
        "projector method, the matrices to convert ('want'; Q, sign, lnabsdet are always returned) and the G fill mode.");
  m.def("dimension_jumps", &dimension_jumps, py::arg("params"), py::arg("E_min"), py::arg("E_max"),
        "Sorted E_cm strictly inside (E_min, E_max) where the basis dimension changes (prompts/002 section 2.2).");
  m.def("compare_kinematics", &compare_kinematics, py::arg("params"),
        "Test hook: oracle make_kinematics vs fast_kinematics, compared bit for bit (momenta, n, ell, m, order).");
  m.def("oracle_point", &oracle_point, py::arg("params"),
        "Oracle (existing CPU path, unchanged physics) at one energy on one core; see qc3.oracle.");
}
