// fast_point.hpp : the cache-free fast path (prompts/001 step 6 rungs 1-4; prompts/002 7.4 K3 pair kernel, 7.5 G pair kernel,
// 7.6 C++ tail).
//
// Same physics, same basis (threshold builder fast_kinematics, bit-identical to qc3_oracle::make_kinematics) and the oracle's
// K2^{-1}; K_df,3 and G use the pair kernels of fast_k3.hpp / fast_g.hpp (same formulas as K3_functions_2plus1.hpp /
// G_functions_v2.h, restructured per momentum pair); the tail (F3 from the block structure of F, Q by one LU solve with N_proj
// right-hand sides, (sign, ln|det|) by LU) is in this file (Eigen, one thread).  The expensive F lattice sum is reorganised:
//   rung 1  one sum pass per spectator momentum for ALL (l'm'; lm) pairs (the oracle repeats the 41^3 cube 16x);
//   rung 2  ellipsoid-within-cube truncation  r^2 <= R^2  with R^2 chosen so that every omitted term is
//           <= trunc_eps  (per-term bound below); NOT a sphere in n;
//   rung 3  real arithmetic (F, G, K2inv, F3, F3^{-1} are real: harmonics are the real Q0 ones, x^2 real);
//   rung 5  (prompts/002 7.7, fast_fclass.hpp, default; Options::f_classes=false gives the direct sum below) class moment tables
//           per canonical nP and one shared sum per (canonical nP, |p|^2); the direct term sum is kept as the cross-check
//   rung 4  energy-independent pieces cached in memory per process (projector V_sel per basis signature = the
//           existing get_projector_cached_v30q map), symmetric block fills.
//
// F-block for spectator momentum p (rows/cols = its (l,m) list), oracle formula F2_ang_mom/I_sum_ang_mom:
//   F_ab = A * ( (2pi/L)^(l_a+l_b) * S_ab  -  delta_ab * I_l ),   A = H/(16 pi^2 L^4 w_p (E-w_p)),
//   S_ab = sum_{n in [-N,N]^3} Y_a(r) Y_b(r) exp(alpha (x^2-r^2)) / (x^2-r^2),
//   r = n_perp + (n_par - xi nP)/gamma  (r = n if nP = 0),  nP = (P-p) L/2pi.
// In the oracle's coordinates r^2 = d^T M d, d = n - xi*nP, M = 1 + (1/gamma^2 - 1) nhat nhat^T.
//
// Truncation bound: |Y_a Y_b| <= 3 r^2 for l <= 1 (Y_1 = sqrt3 r_i), so every term with r^2 > R^2 satisfies
//   |term| <= 3 r^2 e^{-alpha (r^2-x^2)}/(r^2-x^2) <= 3 R^2 e^{-alpha (R^2-x^2)}/(R^2-x^2)   (decreasing for r^2 > x^2 + 1/alpha),
// and R^2 is the smallest value (fixed-point iteration) for which the right side is <= trunc_eps.
#pragma once

#include "oracle_point.hpp"
#include "fast_projector.hpp"
#include "fast_basis.hpp"
#include "fast_k3.hpp"
#include "fast_g.hpp"
#include "fast_fclass.hpp"

namespace qc3_fast {

using qc3_oracle::Mat;
using MatR = Eigen::MatrixXd;
using namespace k3df_fit_v32f;

struct Options {
    double trunc_eps = 1.0e-16;
    std::string projector = "orbit";   // "orbit" (prompts/002 2.1, default) | "eigen" (existing builder, cross-check)
    // which matrices are converted into the Output (prompts/002 2.4: conversions on demand).  The full inverse F3^{-1} is
    // formed ONLY if want_F3inv (Q never needs it).  The K3 basis matrices and K3 are returned iff want_K3.
    bool want_F = true, want_G = true, want_K2inv = true, want_F3 = true, want_F3inv = true, want_Vsel = true, want_K3 = true;
    bool f_classes = true;             // F lattice sum via class moment tables + shared (nP orbit, |p|^2) sums (prompts/002 2.7); false = direct term sum (cross-check)
    bool g_symmetric = false;          // G pair kernel: every ordered pair (default: bit-faithful to the oracle, +0.2 ms at N=160) or symmetric fill (plan 2.3)
};

struct Stats {
    long long lattice_terms = 0;   // summed terms over all momenta (channel 1 + 2)
    long long lattice_cube = 0;    // terms the oracle's cube would have summed per (l'm';lm) entry
    int projector_cache_hit = 0;
    long long f_classes = 0;       // class terms evaluated (class path)
    long long f_sums = 0;          // distinct (canonical nP, |p|^2) moment sums computed
    long long f_blocks = 0;        // spectator-momentum blocks with a non-zero cutoff
};

// smallest R^2 such that 3 R^2 exp(-alpha (R^2-x2))/(R^2-x2) <= eps   (x2 = Re x^2 may be negative)
inline double truncation_R2(double x2, double alpha, double eps) {
    double R2 = std::max(x2, 0.0) + std::log(1.0 / eps) / alpha;
    for (int it = 0; it < 60; ++it) {
        const double d = R2 - x2;
        const double next = x2 + (std::log(3.0 * R2 / d) + std::log(1.0 / eps)) / alpha;
        if (std::abs(next - R2) < 1e-12 * R2) { R2 = next; break; }
        R2 = next;
    }
    return R2 * (1.0 + 1e-12);
}

// Real harmonics l <= 1 (spherical_harmonics of the oracle: Y00 = 1, Y1,-1 = sqrt3 y, Y10 = sqrt3 z, Y11 = sqrt3 x).
inline double Yreal(int ell, int m, double rx, double ry, double rz) {
    static const double s3 = std::sqrt(3.0);
    if (ell == 0) return 1.0;
    if (ell == 1) return s3 * (m == -1 ? ry : (m == 0 ? rz : rx));
    throw std::runtime_error("fast path: only l <= 1 is implemented");
}

// F matrix of one channel (spectator mass mi, pair masses mj, mk). Returns the N_c x N_c real block.
inline MatR F_channel(const std::vector<std::vector<comp>>& cfg, const std::vector<std::vector<int>>& ncfg, comp En,
                      const std::vector<comp>& total_P, double mi, double mj, double mk, double L, double alpha,
                      double epsilon_h, int nmax, bool Q0norm, double trunc_eps, Stats& st, bool use_classes = true) {
    const int n = int(cfg[0].size());
    MatR F = MatR::Zero(n, n);
    if (!Q0norm) throw std::runtime_error("fast path: only the Q0 normalisation is implemented");
    const comp pi = std::acos(-1.0);
    const comp Lby2pi = comp(L) / (comp(2.0) * pi);
    const comp twopibyL = comp(2.0) * pi / comp(L);
    const double tpl = twopibyL.real();
    struct FCacheEntry { int c[3]; int p2; double x2, gamma, xii; FMoments m; };
    std::vector<FCacheEntry> fcache;
    fcache.reserve(64);   // pointers into it must stay valid while the block is assembled

    int i0 = 0;
    while (i0 < n) {
        // rows of one spectator momentum are contiguous (momentum-first basis ordering)
        int i1 = i0 + 1;
        while (i1 < n && ncfg[0][std::size_t(i1)] == ncfg[0][std::size_t(i0)] && ncfg[1][std::size_t(i1)] == ncfg[1][std::size_t(i0)] &&
               ncfg[2][std::size_t(i1)] == ncfg[2][std::size_t(i0)]) ++i1;
        const int nlm = i1 - i0;
        if (nlm > 4) throw std::runtime_error("fast path: more than 4 (l,m) per momentum");
        int ell[4], mm[4];
        for (int a = 0; a < nlm; ++a) {
            ell[a] = int(cfg[3][std::size_t(i0 + a)].real());
            mm[a] = int(cfg[4][std::size_t(i0 + a)].real());
        }
        std::vector<comp> p = {cfg[0][std::size_t(i0)], cfg[1][std::size_t(i0)], cfg[2][std::size_t(i0)]};
        const comp spec_p = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        const comp sigp = sigma_pvec_based(En, p, mi, total_P);
        const comp cutoff = cutoff_function_1(sigp, mj, mk, epsilon_h);
        const comp omega_p = omega_func(spec_p, mi);
        const comp A = cutoff / (comp(16.0) * pi * pi * comp(L) * comp(L) * comp(L) * comp(L) * omega_p * (En - omega_p));
        if (std::abs(A) == 0.0) { i0 = i1; continue; }   // oracle returns exactly 0 for the whole block

        const comp gamma_c = (En - omega_func(spec_p, mi)) / std::sqrt(sigp);
        const comp x_c = std::sqrt(q2psq_star(sigp, mj, mk)) * Lby2pi;
        const comp xii_c = comp(0.5) * (comp(1.0) + comp(mj * mj - mk * mk) / sigp);
        const double gamma = gamma_c.real();
        const double x2 = (x_c * x_c).real();
        const double xibygamma = (xii_c / gamma_c).real();
        const double xii = xii_c.real();
        const double nPx = ((total_P[0] - p[0]) * Lby2pi).real();
        const double nPy = ((total_P[1] - p[1]) * Lby2pi).real();
        const double nPz = ((total_P[2] - p[2]) * Lby2pi).real();
        const double nPsq = nPx * nPx + nPy * nPy + nPz * nPz;
        const bool boosted = (nPsq != 0.0);
        const double invg_m1 = 1.0 / gamma - 1.0;

        const double R2 = truncation_R2(x2, alpha, trunc_eps);
        // quadratic form r^2 = d^T M d, d = n - c0, M = 1 + k nhat nhat^T, k = 1/gamma^2 - 1; centre c0 = xii * nP
        double nh[3] = {0, 0, 0}, c0[3] = {0, 0, 0};
        double kq = 0.0;
        if (boosted) {
            const double nn = std::sqrt(nPsq);
            nh[0] = nPx / nn; nh[1] = nPy / nn; nh[2] = nPz / nn;
            c0[0] = xii * nPx; c0[1] = xii * nPy; c0[2] = xii * nPz;
            kq = 1.0 / (gamma * gamma) - 1.0;
        }
        double M[3][3];
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) M[a][b] = (a == b ? 1.0 : 0.0) + kq * nh[a] * nh[b];

        // moments S_ab accumulated in long double (Neumaier-grade accuracy without a stored term vector)
        long double S[4][4] = {{0}};
        long long count = 0;
        ++st.f_blocks;
        if (use_classes) {
            // ---- class path (fast_fclass.hpp): canonical nP, shared sums, rotation back ----
            int nPi[3] = {int(std::lround(nPx)), int(std::lround(nPy)), int(std::lround(nPz))};
            if (std::abs(nPx - nPi[0]) > 1e-6 || std::abs(nPy - nPi[1]) > 1e-6 || std::abs(nPz - nPi[2]) > 1e-6)
                throw std::runtime_error("fast path: nP is not an integer vector");
            int idx[3] = {0, 1, 2};
            std::stable_sort(idx, idx + 3, [&](int a, int b) { return std::abs(nPi[a]) > std::abs(nPi[b]); });
            int cvec[3], sgn[3];
            for (int j = 0; j < 3; ++j) { cvec[j] = std::abs(nPi[idx[j]]); sgn[j] = (nPi[idx[j]] < 0) ? -1 : 1; }
            const int p2 = ncfg[0][std::size_t(i0)] * ncfg[0][std::size_t(i0)] + ncfg[1][std::size_t(i0)] * ncfg[1][std::size_t(i0)] +
                           ncfg[2][std::size_t(i0)] * ncfg[2][std::size_t(i0)];
            const FMoments* mom = nullptr;
            for (const FCacheEntry& e : fcache)
                if (e.c[0] == cvec[0] && e.c[1] == cvec[1] && e.c[2] == cvec[2] && e.p2 == p2 &&
                    std::abs(e.x2 - x2) <= 1e-10 * (1.0 + std::abs(x2)) && std::abs(e.gamma - gamma) <= 1e-10 * gamma &&
                    std::abs(e.xii - xii) <= 1e-10 * (1.0 + std::abs(xii))) { mom = &e.m; break; }
            if (!mom) {
                FCacheEntry e;
                e.c[0] = cvec[0]; e.c[1] = cvec[1]; e.c[2] = cvec[2]; e.p2 = p2; e.x2 = x2; e.gamma = gamma; e.xii = xii;
                const std::shared_ptr<const FTable> tab = ftable_cached(cvec, nmax, alpha);
                e.m = class_moments(*tab, x2, gamma, xii, R2);
                fcache.push_back(e);
                mom = &fcache.back().m;
                ++st.f_sums;
                st.f_classes += mom->classes;
            }
            // rotate the canonical-frame moments back: r_act[idx[j]] = sgn[j] r_can[j]
            double T1[3], T2[3][3];
            for (int j = 0; j < 3; ++j) T1[idx[j]] = sgn[j] * mom->T1[j];
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) T2[idx[j]][idx[k]] = sgn[j] * sgn[k] * mom->T2[j][k];
            double u[4], v[4][3];
            for (int a = 0; a < nlm; ++a) {
                u[a] = 0.0; v[a][0] = v[a][1] = v[a][2] = 0.0;
                if (ell[a] == 0) u[a] = 1.0;
                else if (ell[a] == 1 && mm[a] >= -1 && mm[a] <= 1) v[a][mm[a] == -1 ? 1 : (mm[a] == 0 ? 2 : 0)] = std::sqrt(3.0);
                else throw std::runtime_error("fast path: only l <= 1 is implemented");
            }
            for (int a = 0; a < nlm; ++a)
                for (int b = a; b < nlm; ++b) {
                    double s = u[a] * u[b] * mom->T0;
                    for (int i = 0; i < 3; ++i) s += (u[a] * v[b][i] + u[b] * v[a][i]) * T1[i];
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j) s += v[a][i] * T2[i][j] * v[b][j];
                    S[a][b] = s;
                }
            count = mom->terms;
        } else {
        const int lo = -nmax, hi = nmax;
        // bounding box in n: half-width R sqrt(1 + (gamma^2-1) nh_i^2) = R sqrt(M^{-1}_ii)
        double hb[3];
        for (int a = 0; a < 3; ++a) hb[a] = std::sqrt(R2) * std::sqrt(1.0 + (gamma * gamma - 1.0) * nh[a] * nh[a]);
        if (!boosted) for (int a = 0; a < 3; ++a) hb[a] = std::sqrt(R2);
        const int xlo = std::max(lo, int(std::ceil(c0[0] - hb[0])) - 1), xhi = std::min(hi, int(std::floor(c0[0] + hb[0])) + 1);
        const int ylo = std::max(lo, int(std::ceil(c0[1] - hb[1])) - 1), yhi = std::min(hi, int(std::floor(c0[1] + hb[1])) + 1);
        for (int ix = xlo; ix <= xhi; ++ix) {
            const double dx = ix - c0[0];
            for (int iy = ylo; iy <= yhi; ++iy) {
                const double dy = iy - c0[1];
                // r^2(dz) = m_zz dz^2 + 2 beta dz + g0 ;  solve <= R2 for dz
                const double beta = M[2][0] * dx + M[2][1] * dy;
                const double g0 = M[0][0] * dx * dx + 2.0 * M[0][1] * dx * dy + M[1][1] * dy * dy;
                const double mzz = M[2][2];
                const double disc = beta * beta - mzz * (g0 - R2);
                if (disc < 0.0) continue;
                const double sq = std::sqrt(disc);
                const double dz_lo = (-beta - sq) / mzz, dz_hi = (-beta + sq) / mzz;
                const int zlo = std::max(lo, int(std::ceil(c0[2] + dz_lo)) - 1), zhi = std::min(hi, int(std::floor(c0[2] + dz_hi)) + 1);
                for (int iz = zlo; iz <= zhi; ++iz) {
                    // r exactly as in I_sum_ang_mom
                    double rx = ix, ry = iy, rz = iz;
                    if (boosted) {
                        const double dot = ix * nPx + iy * nPy + iz * nPz;
                        const double prod1 = (dot / nPsq) * invg_m1 - xibygamma;
                        rx = ix + nPx * prod1; ry = iy + nPy * prod1; rz = iz + nPz * prod1;
                    }
                    const double rsq = rx * rx + ry * ry + rz * rz;
                    if (rsq > R2) continue;
                    const double dd = x2 - rsq;
                    const double w = std::exp(alpha * dd) / dd;
                    double Y[4];
                    for (int a = 0; a < nlm; ++a) Y[a] = Yreal(ell[a], mm[a], rx, ry, rz);
                    for (int a = 0; a < nlm; ++a) {
                        const double wa = w * Y[a];
                        for (int b = a; b < nlm; ++b) S[a][b] += wa * Y[b];
                    }
                    ++count;
                }
            }
        }
        }
        st.lattice_terms += count;
        st.lattice_cube += (long long)(2 * nmax + 1) * (2 * nmax + 1) * (2 * nmax + 1);

        for (int a = 0; a < nlm; ++a) {
            for (int b = a; b < nlm; ++b) {
                const double Bab = double(S[a][b]) * std::pow(tpl, double(ell[a] + ell[b]));
                double Cab = 0.0;
                if (a == b)
                    Cab = I_int_ang_mom(En, sigp, spec_p, std::sqrt(total_P[0] * total_P[0] + total_P[1] * total_P[1] + total_P[2] * total_P[2]),
                                        ell[a], mm[a], ell[b], mm[b], alpha, mi, mj, mk, L, Q0norm).real();
                const double val = (A * comp(Bab - Cab)).real();
                F(i0 + a, i0 + b) = val;
                F(i0 + b, i0 + a) = val;
            }
        }
        i0 = i1;
    }
    return F;
}

inline Mat to_complex(const MatR& M) { return M.cast<comp>(); }

inline double max_abs_imag(const Mat& M) { return M.imag().cwiseAbs().maxCoeff(); }

// (sign, ln|det|) from a partial-pivot LU (invariant 2), same expressions as the oracle's tail.  Singular or non-finite -> domain_error.
template <class M>
inline void slogdet_lu(const M& Q, int& sign, double& sign_imag, double& lnabs) {
    if (!Q.allFinite()) throw std::domain_error("fast: Q contains non-finite entries");
    Eigen::PartialPivLU<M> luQ(Q);
    const M& LU = luQ.matrixLU();
    lnabs = 0.0;
    comp phase(double(luQ.permutationP().determinant()), 0.0);
    for (int k = 0; k < LU.rows(); ++k) {
        const double a = std::abs(LU(k, k));
        if (!(a > 0.0) || !std::isfinite(a)) throw std::domain_error("fast: Q is singular (zero pivot in LU)");
        lnabs += std::log(a);
        phase *= comp(LU(k, k)) / a;
    }
    sign = (phase.real() >= 0.0) ? 1 : -1;
    sign_imag = std::abs(phase.imag());
}

inline qc3_oracle::Output evaluate(const qc3_oracle::Input& in, const Options& opt) {
    using namespace qc3_oracle;
    SingleThread single_thread_guard;
    Output out;
    Clock total, clk;
    out.Ecm = in.Ecm;

    const Kin kin = qc3_fast::fast_kinematics(in);   // threshold basis builder (prompts/002 2.2); the oracle keeps make_kinematics
    const PhysicsParams& par = kin.par;
    out.nnP = kin.spec.nnP;
    out.E_lab = kin.En_c.real();
    out.N1 = kin.N1; out.N2 = kin.N2; out.N = kin.N;
    out.basis = kin.basis;
    const int N = kin.N, N1 = kin.N1, N2 = kin.N2;
    const std::vector<Group> groups = make_groups(kin);
    out.timings["basis"] = clk.lap();

    // ---- F (rungs 1-3) ----
    Stats st;
    MatR Fr = MatR::Zero(N, N);
    Fr.topLeftCorner(N1, N1) = F_channel(kin.plm_config, kin.np_config, kin.En_c, kin.total_P, par.atmK, par.atmK, par.atmpi, kin.L,
                                         par.alpha, par.epsilon_h, int(par.max_shell_num), par.Q0norm, opt.trunc_eps, st, opt.f_classes);
    Fr.bottomRightCorner(N2, N2) = F_channel(kin.klm_config, kin.nk_config, kin.En_c, kin.total_P, par.atmpi, par.atmK, par.atmK, kin.L,
                                             par.alpha, par.epsilon_h, int(par.max_shell_num), par.Q0norm, opt.trunc_eps, st, opt.f_classes);
    out.timings["F"] = clk.lap();

    // ---- K2inv (honours k2_pwave_form): oracle function, unchanged ----
    auto plm = kin.plm_config; auto klm = kin.klm_config;
    Mat K2inv(N, N);
    K2inv_2plus1_by_form(in.k2_pwave_form, K2inv, par.eta_1, par.eta_2, par.scatter_params_1, par.scatter_params_2, kin.En_c, plm, klm, kin.total_P, par.atmK, par.atmpi, par.epsilon_h, kin.L);
    out.timings["K2inv"] = clk.lap();
    {
        const double sc = max_abs_imag(K2inv);
        const double scr = std::max(K2inv.real().cwiseAbs().maxCoeff(), 1e-300);
        if (sc > 1e-14 * scr) throw std::runtime_error("fast: K2inv has a non-negligible imaginary part (real arithmetic invalid)");
    }
    const MatR K2r = K2inv.real();

    // ---- G: pair kernel (prompts/002 2.3, fast_g.hpp), real ----
    const MatR Gr = g_matrix(kin, groups, opt.g_symmetric);
    out.timings["G"] = clk.lap();

    // ---- F3 = F/3 - F H^{-1} F in real arithmetic; F is block diagonal (one block per spectator momentum) ----
    const MatR H = K2r + Fr + Gr;
    Eigen::PartialPivLU<MatR> luH(H);
    const MatR X1 = luH.solve(Fr);
    MatR FX1 = MatR::Zero(N, N);
    for (const Group& g : groups) FX1.middleRows(g.row0, g.nlm).noalias() = Fr.block(g.row0, g.row0, g.nlm, g.nlm) * X1.middleRows(g.row0, g.nlm);
    const MatR F3 = Fr / 3.0 - FX1;
    out.timings["F3"] = clk.lap();

    // ---- projector: fast-owned in-memory cache with the complete key (the legacy map is not used) ----
    bool phit = false;
    FastProjector fproj = fast_projector_cached(kin, opt.projector, phit);
    st.projector_cache_hit = phit ? 1 : 0;
    const Mat& Vsel = fproj.Vsel;
    out.Nproj = fproj.Nproj;
    out.timings["projector"] = clk.lap();
    if (out.Nproj <= 0) throw std::runtime_error("fast: ZERO_PROJECTED_DIM");

    // ---- scaled partial-pivot LU of F3 (the inverse is what gets projected: invariant 1) ----
    const double f3_scale = F3.cwiseAbs().maxCoeff();
    if (!(f3_scale > 0.0) || !std::isfinite(f3_scale)) throw std::runtime_error("fast: BAD_F3_SCALE");
    const MatR F3s = F3 / f3_scale;
    Eigen::PartialPivLU<MatR> luF3(F3s);
    MatR F3invr;
    if (opt.want_F3inv) {          // the full inverse is formed ONLY when requested
        F3invr = luF3.solve(MatR::Identity(N, N));
        F3invr /= f3_scale;
    }
    out.timings["F3inv"] = clk.lap();
    if (!Fr.allFinite() || !F3.allFinite() || !Vsel.allFinite() || (opt.want_F3inv && !F3invr.allFinite())) throw std::runtime_error("fast: NONFINITE_COMPONENT");

    // ---- K3: pair kernel for the four basis matrices (prompts/002 2.3, fast_k3.hpp), K3(c) = sum_i c_i B_i ----
    MatR K3r = MatR::Zero(N, N);
    {
        bool need[4];
        for (int i = 0; i < 4; ++i) need[i] = opt.want_K3 || in.c[std::size_t(i)] != 0.0;
        if (need[0] || need[1] || need[2] || need[3]) {
            std::array<MatR, 4> B;
            k3_basis_matrices(kin, groups, need, B);
            for (int i = 0; i < 4; ++i)
                if (in.c[std::size_t(i)] != 0.0) K3r += in.c[std::size_t(i)] * B[std::size_t(i)];
            if (opt.want_K3) {
                for (int i = 0; i < 4; ++i) out.K3_basis[std::size_t(i)] = to_complex(B[std::size_t(i)]);
                out.has_k3_basis = true;
            }
        }
    }
    out.timings["K3"] = clk.lap();

    // ---- projection of the inverse: Q = L^-6 V^dag (X/s + K3 V) with X = (F3/s)^{-1} V  (N_proj right-hand sides) ----
    // identical to L^-6 V^dag (F3^{-1} + K3) V (invariant 1: project the inverse, scale the matrix, no N_proj factor)
    const double L6 = std::pow(kin.L, 6);
    const bool v_real = (Vsel.imag().cwiseAbs().maxCoeff() == 0.0);
    Mat Q;
    MatR Qr;
    if (v_real) {
        const MatR Vr = Vsel.real();
        const MatR X = luF3.solve(Vr);
        const MatR W = X / f3_scale + K3r * Vr;
        Qr = (Vr.transpose() * W) / L6;
        Q = to_complex(Qr);
    } else {
        MatR rhs(N, 2 * out.Nproj);
        rhs.leftCols(out.Nproj) = Vsel.real();
        rhs.rightCols(out.Nproj) = Vsel.imag();
        const MatR XX = luF3.solve(rhs);
        Mat X(N, out.Nproj);
        X.real() = XX.leftCols(out.Nproj);
        X.imag() = XX.rightCols(out.Nproj);
        const Mat W = X / comp(f3_scale, 0.0) + K3r.cast<comp>() * Vsel;
        Q = (Vsel.adjoint() * W) / comp(L6, 0.0);
    }
    out.timings["projection"] = clk.lap();

    // ---- (sign, ln|det Q|) by partial-pivot LU (invariant 2) ----
    if (v_real) slogdet_lu(Qr, out.sign, out.sign_imag, out.lnabsdet);
    else slogdet_lu(Q, out.sign, out.sign_imag, out.lnabsdet);
    out.timings["lu"] = clk.lap();

    // ---- conversions on demand ----
    if (opt.want_F) out.F = to_complex(Fr);
    if (opt.want_G) out.G = to_complex(Gr);
    if (opt.want_K2inv) out.K2inv = std::move(K2inv);
    if (opt.want_F3) out.F3 = to_complex(F3);
    if (opt.want_F3inv) out.F3inv = to_complex(F3invr);
    if (opt.want_Vsel) out.Vsel = Vsel;
    if (opt.want_K3) out.K3 = to_complex(K3r);
    out.Q = std::move(Q);
    out.timings["total"] = total.lap();
    out.stats["lattice_terms"] = double(st.lattice_terms);
    out.stats["lattice_cube_per_entry"] = double(st.lattice_cube);
    out.stats["f_class_terms"] = double(st.f_classes);
    out.stats["f_sums"] = double(st.f_sums);
    out.stats["f_blocks"] = double(st.f_blocks);
    out.stats["projector_cache_hit"] = double(st.projector_cache_hit);
    return out;
}

}  // namespace qc3_fast
