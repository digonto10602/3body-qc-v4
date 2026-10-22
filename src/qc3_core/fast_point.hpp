// fast_point.hpp : the cache-free fast path v1 (prompts/001 step 6, rungs 1-4).
//
// Same physics, same basis (threshold builder fast_kinematics, bit-identical to qc3_oracle::make_kinematics) and same G, K2^{-1}, K3 functions as the oracle;
// the expensive F lattice sum is reorganised:
//   rung 1  one sum pass per spectator momentum for ALL (l'm'; lm) pairs (the oracle repeats the 41^3 cube 16x);
//   rung 2  ellipsoid-within-cube truncation  r^2 <= R^2  with R^2 chosen so that every omitted term is
//           <= trunc_eps  (per-term bound below); NOT a sphere in n;
//   rung 3  real arithmetic (F, G, K2inv, F3, F3^{-1} are real: harmonics are the real Q0 ones, x^2 real);
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

namespace qc3_fast {

using qc3_oracle::Mat;
using MatR = Eigen::MatrixXd;
using namespace k3df_fit_v32f;

struct Options {
    double trunc_eps = 1.0e-16;
    std::string projector = "orbit";   // "orbit" (prompts/002 2.1, default) | "eigen" (existing builder, cross-check)
};

struct Stats {
    long long lattice_terms = 0;   // summed terms over all momenta (channel 1 + 2)
    long long lattice_cube = 0;    // terms the oracle's cube would have summed per (l'm';lm) entry
    int projector_cache_hit = 0;
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
                      double epsilon_h, int nmax, bool Q0norm, double trunc_eps, Stats& st) {
    const int n = int(cfg[0].size());
    MatR F = MatR::Zero(n, n);
    if (!Q0norm) throw std::runtime_error("fast path: only the Q0 normalisation is implemented");
    const comp pi = std::acos(-1.0);
    const comp Lby2pi = comp(L) / (comp(2.0) * pi);
    const comp twopibyL = comp(2.0) * pi / comp(L);
    const double tpl = twopibyL.real();

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
        st.lattice_terms += count;
        st.lattice_cube += (long long)(hi - lo + 1) * (hi - lo + 1) * (hi - lo + 1);

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
    out.timings["basis"] = clk.lap();

    // ---- F (rungs 1-3) ----
    Stats st;
    MatR Fr = MatR::Zero(N, N);
    Fr.topLeftCorner(N1, N1) = F_channel(kin.plm_config, kin.np_config, kin.En_c, kin.total_P, par.atmK, par.atmK, par.atmpi, kin.L,
                                         par.alpha, par.epsilon_h, int(par.max_shell_num), par.Q0norm, opt.trunc_eps, st);
    Fr.bottomRightCorner(N2, N2) = F_channel(kin.klm_config, kin.nk_config, kin.En_c, kin.total_P, par.atmpi, par.atmK, par.atmK, kin.L,
                                             par.alpha, par.epsilon_h, int(par.max_shell_num), par.Q0norm, opt.trunc_eps, st);
    out.timings["F"] = clk.lap();

    // ---- K2inv (honours k2_pwave_form), G: oracle functions, unchanged ----
    auto plm = kin.plm_config; auto klm = kin.klm_config;
    Mat K2inv(N, N), G(N, N);
    K2inv_2plus1_by_form(in.k2_pwave_form, K2inv, par.eta_1, par.eta_2, par.scatter_params_1, par.scatter_params_2, kin.En_c, plm, klm, kin.total_P, par.atmK, par.atmpi, par.epsilon_h, kin.L);
    out.timings["K2inv"] = clk.lap();
    G_2plus1_mat(G, kin.En_c, plm, klm, kin.total_P, par.atmK, par.atmpi, kin.L, par.alpha, par.epsilon_h, par.max_shell_num, par.Q0norm);
    out.timings["G"] = clk.lap();

    // ---- F3 = F/3 - F H^{-1} F in real arithmetic ----
    {
        const double sc = std::max({max_abs_imag(K2inv), max_abs_imag(G)});
        const double scr = std::max({K2inv.real().cwiseAbs().maxCoeff(), G.real().cwiseAbs().maxCoeff(), 1e-300});
        if (sc > 1e-14 * scr) throw std::runtime_error("fast: K2inv or G has a non-negligible imaginary part (real arithmetic invalid)");
    }
    const MatR K2r = K2inv.real(), Gr = G.real();
    const MatR H = K2r + Fr + Gr;
    Eigen::PartialPivLU<MatR> luH(H);
    const MatR X1 = luH.solve(Fr);
    MatR F3 = Fr / 3.0 - Fr * X1;
    out.timings["F3"] = clk.lap();

    // ---- projector: fast-owned in-memory cache with the complete key (the legacy map is not used) ----
    bool phit = false;
    FastProjector fproj = fast_projector_cached(kin, opt.projector, phit);
    st.projector_cache_hit = phit ? 1 : 0;
    Mat Vsel = fproj.Vsel;
    out.Nproj = fproj.Nproj;
    out.timings["projector"] = clk.lap();
    if (out.Nproj <= 0) throw std::runtime_error("fast: ZERO_PROJECTED_DIM");

    // ---- F3^{-1}: scaled partial-pivot LU, full inverse (the inverse is what gets projected: invariant 1) ----
    const double f3_scale = F3.cwiseAbs().maxCoeff();
    if (!(f3_scale > 0.0) || !std::isfinite(f3_scale)) throw std::runtime_error("fast: BAD_F3_SCALE");
    const MatR F3s = F3 / f3_scale;
    Eigen::PartialPivLU<MatR> luF3(F3s);
    MatR F3invr = luF3.solve(MatR::Identity(N, N));
    F3invr /= f3_scale;
    out.timings["F3inv"] = clk.lap();
    if (!Fr.allFinite() || !F3.allFinite() || !F3invr.allFinite() || !Vsel.allFinite()) throw std::runtime_error("fast: NONFINITE_COMPONENT");

    // ---- K3 ----
    const bool any_c = (in.c[0] != 0.0 || in.c[1] != 0.0 || in.c[2] != 0.0 || in.c[3] != 0.0);
    Mat K3 = Mat::Zero(N, N);
    if (any_c || in.compute_k3) {
        compute_k3(kin, in.c, out.K3_basis, K3);
        out.has_k3_basis = true;
    }
    out.timings["K3"] = clk.lap();

    out.F = to_complex(Fr); out.G = std::move(G); out.K2inv = std::move(K2inv);
    out.F3 = to_complex(F3); out.F3inv = to_complex(F3invr); out.Vsel = std::move(Vsel); out.K3 = std::move(K3);
    out.timings["total"] = total.lap();
    out.stats["lattice_terms"] = double(st.lattice_terms);
    out.stats["lattice_cube_per_entry"] = double(st.lattice_cube);
    out.stats["projector_cache_hit"] = double(st.projector_cache_hit);
    return out;
}

}  // namespace qc3_fast
