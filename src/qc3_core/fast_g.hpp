// fast_g.hpp : pair kernel for G (prompts/002 section 2.3, sub-step 7.5).
//
// Physics: EXACTLY G_ij_lm / G_2plus1_mat of src/core/G_functions_v2.h (not edited):
//   G_{p l'm' i; k l m j} = Y_{l'm'}(k*) Y_{lm}(p*) cutoff(sig_p) cutoff(sig_k) / (2 w_p L^3 * denom * 2 w_k L^3),
//   k* = boost(w_k, k, E - w_p, P - p),  p* = boost(w_p, p, E - w_k, P - k),
//   denom = (E - w_p - w_k)^2 - |P - p - k|^2 - m_third^2,
//   flavour factors sqrt2 (-1)^{l of the kaon row/column} for (1,2) and (2,1), and a zero (2,2) block (Q0 normalisation:
//   no (q*)^l factors).
// Restructured: per-momentum scalars (omega, E-omega, P-p, boost parameters, 1/(2 w L^3), sigma and the cutoffs for each
// partner flavour) are computed once per momentum with the oracle's own helpers (comp); per pair, k* and p* use the stored
// boost parameters and the whole (nlm_A x nlm_B) block is filled.
// Symmetry: G(p,k) = G(k,p) mathematically (the same Y's, cutoffs, denominator and flavour factor).  Numerically the oracle's
// matrix is symmetric only to <= 5.4e-14 relative Frobenius (worst case 200_A2 0.34738: the denominator
// (E - w_p) - w_k vs (E - w_k) - w_p cancels near its zero).  Two fill modes:
//   symmetric = true  : pairs A <= B plus the transposed block (the plan's "symmetric fill"); the mirrored half then differs
//                       from the oracle by that asymmetry, and ln|det Q| by up to 1.3e-8 at 100_A2 0.2915561 (tolerance 1e-6).
//   symmetric = false : every ordered pair, same operation order as the oracle: G is bit-identical to G_2plus1_mat (measured
//                       deviation 0.0 at the 12 anchors) for +0.2 ms at N = 160.  This is the default of the fast path.
#pragma once

#include "fast_pairs.hpp"

namespace qc3_fast {

struct GMom {
    int flav = 0;
    double p[3] = {0, 0, 0};
    double w = 0.0;            // omega_func(|p|, m)
    double E2 = 0.0;           // E - w
    double P2[3] = {0, 0, 0};  // P - p
    BoostPar bp;               // boost parameters of (E2, P2)
    double inv2wL3 = 0.0;      // 1/(2 w L^3)
    double cut[3] = {0, 0, 0}; // [v]: cutoff_function_1 with the (mj, mk) the oracle uses when the partner has flavour v
};

inline void g_momentum(const qc3_oracle::Kin& kin, const Group& g, GMom& r) {
    const double M1 = kin.par.atmK, M2 = kin.par.atmpi, L = kin.L, eps_h = kin.par.epsilon_h;
    const comp En = kin.En_c;
    const double mi = (g.flav == 1) ? M1 : M2;
    r.flav = g.flav;
    for (int a = 0; a < 3; ++a) r.p[a] = g.p[a];
    const comp spec = std::sqrt(g.pc[0] * g.pc[0] + g.pc[1] * g.pc[1] + g.pc[2] * g.pc[2]);
    const comp w = omega_func(spec, mi);
    const comp E2 = En - w;
    comp P2[3];
    for (int a = 0; a < 3; ++a) P2[a] = kin.total_P[std::size_t(a)] - g.pc[a];
    const std::vector<comp> pv = {g.pc[0], g.pc[1], g.pc[2]};
    const comp sig = sigma_pvec_based(En, pv, mi, kin.total_P);
    check_real(w, "omega"); check_real(E2, "E-omega"); check_real(sig, "sigma");
    r.w = w.real(); r.E2 = E2.real();
    for (int a = 0; a < 3; ++a) { check_real(P2[a], "P-p"); r.P2[a] = P2[a].real(); }
    r.bp = make_boost_par(P2, E2);
    const comp o = 1.0 / (2.0 * w * L * L * L);      // oneby2omegapLcube
    check_real(o, "1/(2 w L^3)");
    r.inv2wL3 = o.real();
    for (int v = 1; v <= 2; ++v) {
        double mj, mk;     // arguments (mj, mk) of cutoff_function_1(sig, mj, mk, eps_h) in G_ij_lm for this momentum
        if (g.flav == 1) { mj = (v == 1 ? M1 : M2); mk = (v == 1 ? M2 : M1); }   // G11: (m1, m2); G12/G21: (m2, m1)
        else { mj = M1; mk = M1; }
        const comp cu = cutoff_function_1(sig, mj, mk, eps_h);
        check_real(cu, "cutoff");
        r.cut[v] = cu.real();
    }
}

// Block of G for the ORDERED pair (A as row momentum p, B as column momentum k), including the flavour factor.
inline void g_pair_block(const qc3_oracle::Kin& kin, const GMom& A, const GMom& B, const Group& gA, const Group& gB, double Gv[4][4]) {
    const double M1 = kin.par.atmK, M2 = kin.par.atmpi;
    const double E = kin.En_c.real();
    const double* P = nullptr;
    double Pd[3];
    for (int a = 0; a < 3; ++a) Pd[a] = kin.total_P[std::size_t(a)].real();
    P = Pd;
    double kstar[3], pstar[3];
    boost_apply(A.bp, B.w, B.p, kstar);       // k_star = boost(omg_k, k, En - omg_p, P - p)
    boost_apply(B.bp, A.w, A.p, pstar);       // p_star = boost(omg_p, p, En - omg_k, P - k)
    const double mk = (A.flav == 1 && B.flav == 1) ? M2 : M1;   // third particle: G11 -> m2 ; G12, G21 -> m1
    // denom = (En - w_p - w_k)*(En - w_p - w_k) - mom_P_p_k*mom_P_p_k - mk*mk,  mom_P_p_k = sqrt(|P - p - k|^2)
    const double t = E - A.w - B.w;
    const double mom = std::sqrt((P[0] - A.p[0] - B.p[0]) * (P[0] - A.p[0] - B.p[0]) + (P[1] - A.p[1] - B.p[1]) * (P[1] - A.p[1] - B.p[1]) +
                                 (P[2] - A.p[2] - B.p[2]) * (P[2] - A.p[2] - B.p[2]));
    const double denom = t * t - mom * mom - mk * mk;
    const double onebydenom = 1.0 / denom;
    const double c1 = A.cut[B.flav], c2 = B.cut[A.flav];
    static const double sqrt2 = std::sqrt(2.0);
    for (int a = 0; a < gA.nlm; ++a)
        for (int b = 0; b < gB.nlm; ++b) {
            const double Y1 = Yq0(gA.ell[a], gA.m[a], kstar);   // Ylm1 = Y_{l_f m_f}(k*)
            const double Y2 = Yq0(gB.ell[b], gB.m[b], pstar);   // Ylm2 = Y_{l_i m_i}(p*)
            // result = Ylm1*oneby2omegapLcube*cutoff1*cutoff2*onebydenom*oneby2omegakLcube*Ylm2
            double r = Y1 * A.inv2wL3;
            r = r * c1; r = r * c2; r = r * onebydenom; r = r * B.inv2wL3; r = r * Y2;
            if (A.flav == 1 && B.flav == 2) r = sqrt2 * ((gA.ell[a] % 2 != 0) ? -1.0 : 1.0) * r;     // sqrt(2) * (-1)^{l_f} * G
            else if (A.flav == 2 && B.flav == 1) r = sqrt2 * r * ((gB.ell[b] % 2 != 0) ? -1.0 : 1.0);   // sqrt(2) * G * (-1)^{l_i}
            Gv[a][b] = r;
        }
}

// Full N x N real G.  symmetric = true: pairs A <= B (and the transposed block); false: every ordered pair.
inline MatR g_matrix(const qc3_oracle::Kin& kin, const std::vector<Group>& gs, bool symmetric) {
    const int N = kin.N;
    const int ng = int(gs.size());
    std::vector<GMom> mom(static_cast<std::size_t>(ng));
    for (int g = 0; g < ng; ++g) g_momentum(kin, gs[std::size_t(g)], mom[std::size_t(g)]);
    MatR Gm = MatR::Zero(N, N);
    double Gv[4][4];
    for (int ia = 0; ia < ng; ++ia) {
        const Group& gA = gs[std::size_t(ia)];
        for (int ib = (symmetric ? ia : 0); ib < ng; ++ib) {
            const Group& gB = gs[std::size_t(ib)];
            if (gA.flav == 2 && gB.flav == 2) continue;      // zero (2,2) block (Filler0_22)
            g_pair_block(kin, mom[std::size_t(ia)], mom[std::size_t(ib)], gA, gB, Gv);
            for (int a = 0; a < gA.nlm; ++a)
                for (int b = 0; b < gB.nlm; ++b) {
                    if (symmetric && ia == ib && b < a) continue;                  // diagonal block: upper triangle only
                    Gm(gA.row0 + a, gB.row0 + b) = Gv[a][b];
                    if (symmetric) Gm(gB.row0 + b, gA.row0 + a) = Gv[a][b];
                }
        }
    }
    if (!Gm.allFinite()) throw std::runtime_error("fast: nonfinite G (pair kernel)");
    return Gm;
}

}  // namespace qc3_fast
