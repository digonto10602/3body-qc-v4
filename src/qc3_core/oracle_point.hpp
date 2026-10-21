// oracle_point.hpp : the ORACLE of prompts/001 -- the existing CPU path, unchanged, at one energy, one core.
//
// This file only orchestrates calls into the existing physics headers, in exactly the order of
// build_component_entry() (src/cachegen_cpu/generate_cpu_component_cache.cpp, L223-288):
//   config_maker_4_momentum_first (basis) -> F2_2plus1_mat -> K2inv_EREord2_2plus1_mat -> G_2plus1_mat
//   -> H = K2inv + F2 + G -> F3 = F2/3 - F2 H^{-1} F2 -> get_projector_cached_v30q (Vsel)
//   -> scaled partial-pivot LU F3inv -> V^dag F3inv V,
// plus k3_2plus1::K3mat_2plus1 for the four basis matrices (c = e_i), as in precompute_projected_k3_basis().
// No formula is copied.  Final assembly (CLAUDE.md invariants 1-2):
//   Q = (L)^-6 V^dag (F3inv + K3) V,  L = xi*Lbyas;  (sign, ln|det Q|) by partial-pivot LU.
//
// The projector cache of the existing code is a process-global map keyed by the basis signature; it is cleared
// at the start of every call so that one oracle call really is "one energy from scratch" (the cost that the
// 2.8-6.5 s of research 007 refers to).
#pragma once

#define V32F_DISABLE_MINUIT 1
#define V31L_DISABLE_MINUIT 1
#define V31R_DISABLE_MINUIT 1
#define V31S_DISABLE_MINUIT 1
#define V31T_DISABLE_MINUIT 1
#ifndef V32F_HAS_MINUIT2
#define V32F_HAS_MINUIT2 0
#endif

#include "K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp"
#include "k2inv_forms.hpp"

#include <Eigen/Dense>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace qc3_oracle {

using k3df_fit_v32f::FitSettings;
using Mat = Eigen::MatrixXcd;

struct Input {
    std::string irrep;               // "000_A1m", "100_A2", ... (parse_label of the existing code)
    double Lbyas = 20.0, xi = 3.444, atmK = 0.09698, atmpi = 0.06906;
    double eta_1 = 1.0, eta_2 = 0.5, alpha = 0.5, epsilon_H = 0.0;
    int max_shell_num = 20;
    std::vector<int> waves_1 = {0, 1}, waves_2 = {0};
    int parity = -1;
    double eig_tol = 0.05;
    double a0_Kpi = 4.04, a1_Kpi = -43.2, a0_KK = 4.12;
    std::string k2_pwave_form = "coded_v2";   // "coded_v2" (production) | "q0_consistent" (prompts/001 section 2.2)
    double Ecm = 0.0;
    std::array<double, 4> c = {0.0, 0.0, 0.0, 0.0};
    bool compute_k3 = true;          // evaluate the four K3 basis matrices (always done if any c != 0)
};

struct Output {
    double Ecm = 0, E_lab = 0;
    int N1 = 0, N2 = 0, N = 0, Nproj = 0;
    std::array<int, 3> nnP{{0, 0, 0}};
    std::vector<std::array<int, 6>> basis;   // flavour, nx, ny, nz, ell, m
    Mat F, G, K2inv, F3, F3inv, Vsel, K3, Q;
    std::array<Mat, 4> K3_basis;
    bool has_k3_basis = false;
    int sign = 0;
    double sign_imag = 0.0;                  // |imaginary part of the LU phase product|, should be ~1e-16
    double lnabsdet = 0.0;
    std::map<std::string, double> timings;
};

class Clock {
    std::chrono::steady_clock::time_point t0_ = std::chrono::steady_clock::now();
public:
    double lap() {
        const auto t1 = std::chrono::steady_clock::now();
        const double s = std::chrono::duration<double>(t1 - t0_).count();
        t0_ = t1;
        return s;
    }
};

// RAII: force one OpenMP thread and one Eigen thread inside the call, restore afterwards.
struct SingleThread {
#ifdef _OPENMP
    int omp_prev;
#endif
    int eigen_prev;
    SingleThread() {
#ifdef _OPENMP
        omp_prev = omp_get_max_threads();
        omp_set_num_threads(1);
#endif
        eigen_prev = Eigen::nbThreads();
        Eigen::setNbThreads(1);
    }
    ~SingleThread() {
#ifdef _OPENMP
        omp_set_num_threads(omp_prev);
#endif
        Eigen::setNbThreads(eigen_prev);
    }
};

inline Output evaluate(const Input& in) {
    using namespace k3df_fit_v32f;
    SingleThread single_thread_guard;
    Output out;
    Clock total, clk;
    out.Ecm = in.Ecm;

    // ---- parameters: same objects the cache generator builds (FitSettings -> make_base_physics) ----
    FitSettings s;
    s.Lval = in.Lbyas; s.xival = in.xi;
    s.atmpi = in.atmpi; s.atmK = in.atmK;
    s.eta_1 = in.eta_1; s.eta_2 = in.eta_2;
    s.alpha = in.alpha; s.epsilon_h = in.epsilon_H;
    s.max_shell_num = double(in.max_shell_num);
    s.parity = in.parity; s.eig_tol = in.eig_tol;
    s.waves_vec_1 = in.waves_1; s.waves_vec_2 = in.waves_2;
    s.scatter_params_1[0][0] = in.a0_Kpi;
    s.scatter_params_1[1][0] = in.a1_Kpi;
    s.scatter_params_2[0][0] = in.a0_KK;
    s.omp_threads = 1;
    PhysicsParams par = make_base_physics(s);
    par.waves_vec_1 = in.waves_1;
    par.waves_vec_2 = in.waves_2;
    const MomentumIrrepSpec spec = parse_label(in.irrep);
    out.nnP = spec.nnP;

    {   // fresh projector cache: one oracle call = one energy from scratch
        std::lock_guard<std::mutex> lock(g_projector_cache_mutex_v30q);
        g_projector_cache_v30q.clear();
    }

    // ---- kinematics + basis (build_component_entry L226-240) ----
    const comp pi = std::acos(-1.0);
    const double L = par.L();
    const comp twopibyL = comp(2.0, 0.0) * pi / comp(L, 0.0);
    std::vector<comp> total_P(3), nnP_config(3);
    for (int a = 0; a < 3; ++a) {
        total_P[a] = twopibyL * double(spec.nnP[a]);
        nnP_config[a] = comp(spec.nnP[a], 0.0);
    }
    const comp Ecm_c(in.Ecm, 0.0);
    const comp En_c = Ecm_to_E(Ecm_c, total_P);
    out.E_lab = En_c.real();

    std::vector<std::vector<comp>> plm_config(5), klm_config(5);
    std::vector<std::vector<int>> np_config(5), nk_config(5);
    config_maker_4_momentum_first(plm_config, np_config, par.waves_vec_1, En_c, total_P, par.atmK, par.atmK, par.atmpi, L, par.epsilon_h, par.max_shell_num, par.tolerance);
    config_maker_4_momentum_first(klm_config, nk_config, par.waves_vec_2, En_c, total_P, par.atmpi, par.atmK, par.atmK, L, par.epsilon_h, par.max_shell_num, par.tolerance);
    out.N1 = int(plm_config[0].size());
    out.N2 = int(klm_config[0].size());
    out.N = out.N1 + out.N2;
    if (out.N <= 0) throw std::runtime_error("oracle: EMPTY_CONFIG (no basis states at this energy)");
    out.basis.reserve(std::size_t(out.N));
    for (int a = 0; a < out.N1; ++a)
        out.basis.push_back({1, np_config[0][a], np_config[1][a], np_config[2][a], np_config[3][a], np_config[4][a]});
    for (int a = 0; a < out.N2; ++a)
        out.basis.push_back({2, nk_config[0][a], nk_config[1][a], nk_config[2][a], nk_config[3][a], nk_config[4][a]});
    out.timings["basis"] = clk.lap();

    // ---- ingredients ----
    const int N = out.N;
    Mat F2(N, N), G(N, N), K2inv(N, N);
    F2_2plus1_mat(F2, En_c, plm_config, klm_config, total_P, par.atmK, par.atmpi, L, par.alpha, par.epsilon_h, par.max_shell_num, par.Q0norm);
    out.timings["F"] = clk.lap();
    K2inv_2plus1_by_form(in.k2_pwave_form, K2inv, par.eta_1, par.eta_2, par.scatter_params_1, par.scatter_params_2, En_c, plm_config, klm_config, total_P, par.atmK, par.atmpi, par.epsilon_h, L);
    out.timings["K2inv"] = clk.lap();
    G_2plus1_mat(G, En_c, plm_config, klm_config, total_P, par.atmK, par.atmpi, L, par.alpha, par.epsilon_h, par.max_shell_num, par.Q0norm);
    out.timings["G"] = clk.lap();

    const Mat H = K2inv + F2 + G;
    Eigen::PartialPivLU<Mat> luH(H);
    Mat X1 = luH.solve(F2);
    Mat F3 = (F2 / comp(3.0, 0.0)) - F2 * X1;
    out.timings["F3"] = clk.lap();

    CachedProjectorV30q pcache = get_projector_cached_v30q(plm_config, np_config, klm_config, nk_config, spec.irrep, nnP_config, par);
    Mat Vsel = pcache.Vsel;
    out.Nproj = int(Vsel.cols());
    out.timings["projector"] = clk.lap();
    if (out.Nproj <= 0) throw std::runtime_error("oracle: ZERO_PROJECTED_DIM");

    double f3_scale = 0.0;
    for (int r = 0; r < F3.rows(); ++r) for (int c = 0; c < F3.cols(); ++c) {
        const double a = std::abs(F3(r, c));
        if (std::isfinite(a) && a > f3_scale) f3_scale = a;
    }
    if (!(f3_scale > 0.0) || !std::isfinite(f3_scale)) throw std::runtime_error("oracle: BAD_F3_SCALE");
    Mat F3_scaled = F3 / f3_scale;
    Eigen::PartialPivLU<Mat> luF3(F3_scaled);
    Mat F3inv = luF3.solve(Mat::Identity(N, N));
    F3inv /= f3_scale;
    out.timings["F3inv"] = clk.lap();
    if (!K2inv.allFinite() || !F2.allFinite() || !G.allFinite() || !F3.allFinite() || !F3inv.allFinite() || !Vsel.allFinite())
        throw std::runtime_error("oracle: NONFINITE_COMPONENT");

    // ---- K3: four basis matrices (c = e_i), K3(c) = sum c_i B_i ----
    const bool any_c = (in.c[0] != 0.0 || in.c[1] != 0.0 || in.c[2] != 0.0 || in.c[3] != 0.0);
    Mat K3 = Mat::Zero(N, N);
    if (any_c || in.compute_k3) {
        const double M1 = par.atmK, M2 = par.atmpi;
        for (int i = 0; i < 4; ++i) {
            std::vector<comp> Kiso = {comp(i == 0 ? 1.0 : 0.0, 0.0), comp(i == 1 ? 1.0 : 0.0, 0.0)};
            const comp KB = comp(i == 2 ? 1.0 : 0.0, 0.0), KE = comp(i == 3 ? 1.0 : 0.0, 0.0);
            Mat B(N, N);
            k3_2plus1::K3mat_2plus1(B, En_c, plm_config, klm_config, total_P, M1, M2, Kiso, KB, KE, 'n');
            if (!B.allFinite()) throw std::runtime_error("oracle: nonfinite K3 basis matrix");
            out.K3_basis[std::size_t(i)] = B;
            K3 += comp(in.c[std::size_t(i)], 0.0) * B;
        }
        out.has_k3_basis = true;
    }
    out.timings["K3"] = clk.lap();

    // ---- projection: Q = L^-6 V^dag (F3inv + K3) V  (invariant 1: project the inverse, scale the matrix) ----
    Mat Q = (Vsel.adjoint() * (F3inv + K3) * Vsel) / comp(std::pow(L, 6), 0.0);
    out.timings["projection"] = clk.lap();

    // ---- (sign, ln|det Q|) by partial-pivot LU (invariant 2) ----
    {
        Eigen::PartialPivLU<Mat> luQ(Q);
        const Mat& LU = luQ.matrixLU();
        double lnabs = 0.0;
        comp phase(double(luQ.permutationP().determinant()), 0.0);
        for (int k = 0; k < LU.rows(); ++k) {
            const double a = std::abs(LU(k, k));
            lnabs += std::log(a);
            phase *= (a > 0.0) ? LU(k, k) / a : comp(1.0, 0.0);
        }
        out.lnabsdet = lnabs;
        out.sign = (phase.real() >= 0.0) ? 1 : -1;
        out.sign_imag = std::abs(phase.imag());
    }
    out.timings["lu"] = clk.lap();

    out.F = std::move(F2); out.G = std::move(G); out.K2inv = std::move(K2inv);
    out.F3 = std::move(F3); out.F3inv = std::move(F3inv); out.Vsel = std::move(Vsel);
    out.K3 = std::move(K3); out.Q = std::move(Q);
    out.timings["total"] = total.lap();
    return out;
}

}  // namespace qc3_oracle
