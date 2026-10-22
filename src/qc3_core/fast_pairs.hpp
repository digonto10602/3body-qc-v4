// fast_pairs.hpp : shared pieces of the pair kernels for K_df,3 (fast_k3.hpp) and G (fast_g.hpp)
// (prompts/002 section 2.3, sub-steps 7.4 and 7.5).
//
// The basis is ordered momentum-first: the rows of one spectator momentum (flavour f, momentum p) are contiguous and
// carry its (l, m) list.  A "group" is one such momentum.  A pair kernel visits every unordered pair of groups once,
// computes the pair kinematics once, and fills the whole (nlm_A x nlm_B) block of every matrix it produces.
//
// Arithmetic: the oracle evaluates these functions in std::complex<double>, but for real E_lab, real momenta and
// sigma > sigma_min > 0 (guaranteed by the basis threshold) every complex quantity has an exactly zero imaginary part,
// and a product/quotient of two such numbers equals the real operation (the libgcc complex division gives a/c exactly
// for zero imaginary parts).  The per-pair arithmetic is therefore done in double, in the SAME operation order as the
// oracle's formulas (written next to each formula).  The per-momentum scalars (omega, sigma, q*^2, cutoffs, boost
// parameters) are computed ONCE per momentum with the oracle's own scalar helpers (comp arithmetic) and converted.
// Any sqrt of a negative number would yield NaN here (the oracle: an imaginary part); every result matrix is checked
// for finiteness by the caller, so the real-domain assumption cannot fail silently.
#pragma once

#include "oracle_point.hpp"

namespace qc3_fast {

using MatR = Eigen::MatrixXd;

// One spectator momentum of the basis (all its (l, m) rows).
struct Group {
    int flav = 0;          // 1: kaon spectator (channel 1), 2: pion spectator (channel 2)
    int row0 = 0;          // first row in the full N x N basis
    int nlm = 0;
    int ell[4] = {0, 0, 0, 0};
    int m[4] = {0, 0, 0, 0};
    comp pc[3];            // the oracle's complex momentum (imag == 0)
    double p[3] = {0, 0, 0};
};

inline std::vector<Group> make_groups(const qc3_oracle::Kin& kin) {
    std::vector<Group> gs;
    auto add = [&](int flav, int rowoff, const std::vector<std::vector<comp>>& cfg, const std::vector<std::vector<int>>& ncfg) {
        const int n = int(cfg[0].size());
        int i0 = 0;
        while (i0 < n) {
            int i1 = i0 + 1;
            while (i1 < n && ncfg[0][std::size_t(i1)] == ncfg[0][std::size_t(i0)] && ncfg[1][std::size_t(i1)] == ncfg[1][std::size_t(i0)] &&
                   ncfg[2][std::size_t(i1)] == ncfg[2][std::size_t(i0)]) ++i1;
            Group g;
            g.flav = flav; g.row0 = rowoff + i0; g.nlm = i1 - i0;
            if (g.nlm > 4) throw std::runtime_error("fast pair kernels: more than 4 (l,m) per momentum");
            for (int a = 0; a < g.nlm; ++a) {
                g.ell[a] = int(std::llround(cfg[3][std::size_t(i0 + a)].real()));
                g.m[a] = int(std::llround(cfg[4][std::size_t(i0 + a)].real()));
                if (g.ell[a] < 0 || g.ell[a] > 1) throw std::runtime_error("fast pair kernels: only l <= 1 is implemented");
            }
            for (int c = 0; c < 3; ++c) {
                g.pc[c] = cfg[std::size_t(c)][std::size_t(i0)];
                if (g.pc[c].imag() != 0.0) throw std::runtime_error("fast pair kernels: complex momentum");
                g.p[c] = g.pc[c].real();
            }
            gs.push_back(g);
            i0 = i1;
        }
    };
    add(1, 0, kin.plm_config, kin.np_config);
    add(2, kin.N1, kin.klm_config, kin.nk_config);
    return gs;
}

// Parameters of boost(p0, p, E2, P2) (functions.h): P2hat, beta2 = |P2|/E2, gam2 = 1/sqrt(1-beta2^2); exact-zero P2 returns p.
struct BoostPar {
    bool zero = true;
    double hat[3] = {0, 0, 0};
    double beta = 0.0, gam = 1.0;
};

inline BoostPar make_boost_par(const comp P2[3], comp E2) {
    BoostPar b;
    const comp zeroval = {0.0, 0.0};
    if (P2[0] == zeroval && P2[1] == zeroval && P2[2] == zeroval) return b;
    b.zero = false;
    const comp P2norm = std::sqrt(P2[0] * P2[0] + P2[1] * P2[1] + P2[2] * P2[2]);
    const comp beta2 = P2norm / E2;
    const comp gam2 = 1.0 / std::sqrt(1.0 - beta2 * beta2);
    for (int c = 0; c < 3; ++c) b.hat[c] = (P2[c] / P2norm).real();
    b.beta = beta2.real();
    b.gam = gam2.real();
    return b;
}

// out = boost(p0, p, E2, P2) with the parameters of (E2, P2), same operation order as boost() in functions.h.
inline void boost_apply(const BoostPar& b, double p0, const double p[3], double out[3]) {
    if (b.zero) { out[0] = p[0]; out[1] = p[1]; out[2] = p[2]; return; }
    const double dot_p_P2hat = p[0] * b.hat[0] + p[1] * b.hat[1] + p[2] * b.hat[2];
    const double common = (b.gam - 1.0) * dot_p_P2hat - b.gam * b.beta * p0;
    out[0] = p[0] + common * b.hat[0];
    out[1] = p[1] + common * b.hat[1];
    out[2] = p[2] + common * b.hat[2];
}

inline double dot3d(const double a[3], const double b[3]) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

// |p|^l Y_lm(p) for l <= 1 in the Q0 normalisation of spherical_harmonics(): Y00 = 1, Y1,-1 = sqrt3 y, Y10 = sqrt3 z, Y11 = sqrt3 x.
inline double Yq0(int ell, int m, const double v[3]) {
    static const double s3 = std::sqrt(comp(3.0)).real();
    if (ell == 0) return 1.0;
    return s3 * (m == -1 ? v[1] : (m == 0 ? v[2] : v[0]));
}

inline void check_real(comp x, const char* what) {
    if (x.imag() != 0.0) throw std::runtime_error(std::string("fast pair kernels: complex ") + what + " (real arithmetic invalid)");
}

}  // namespace qc3_fast
