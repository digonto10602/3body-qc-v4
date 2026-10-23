// fast_k3.hpp : pair kernel for the four K_df,3 basis matrices (prompts/002 section 2.3, sub-step 7.4).
//
// Physics: EXACTLY the formulas of src/core/K3_functions_2plus1.hpp (K3B_element, K3E_element, K3_element_raw,
// K3mat_2plus1, two_plus_one_block_factor).  That file is not edited (it is the oracle); this kernel restructures the
// evaluation:
//   * everything that depends on ONE momentum (omega, E-omega, P-p, sigma, q*^2 -> pms_0, the boost parameters, beta
//     vectors, gamma factors, p'_1 = boost(...), omega_{p'_1}) is computed once per momentum with the oracle's own
//     scalar helpers (k3_2plus1::omega_vec, sigma_i_cpp, kallentriangle, boost_cpp, operators on Vec3), in comp;
//   * everything that depends on the PAIR (the boosted vectors psk, ksp, P2psk, P2ksp, the vectors V, V', the tensor t,
//     the scalar dot products) is computed once per unordered pair and fills the whole (l'm'; lm) block of ALL
//     FOUR matrices (iso0, iso1, B, E), then the transposed block (the matrices are symmetric by the formulas; measured
//     on the oracle: iso0, iso1 exactly, B <= 1e-17, E <= 8.2e-17 relative Frobenius, prompts/002 section 2.3).
// Re-expression (justified): the per-pair arithmetic is double instead of comp.  All imaginary parts are exactly zero
// for real E_lab, real momenta and sigma > 0 (the basis threshold), so the real operation equals the complex one; the
// operation order of every expression is that of the oracle (the oracle line is quoted next to each formula).
// Verified against K3mat_2plus1 to <= 1e-12 relative (acceptance test_m1_fast_vs_oracle, unit tests/unit/test_fast_k3_g.py).
// Basis pair order: channel 1 (kaon spectator) before channel 2, so a pair (A,B), A <= B, has flavour (1,1), (1,2) or
// (2,2); the (2,1) block is the transpose of the (1,2) block.
#pragma once

#include "fast_pairs.hpp"

namespace qc3_fast {

struct K3Mom {
    int flav = 0;
    double Mi = 0.0;
    double p[3] = {0, 0, 0};
    double om = 0.0;               // omega_vec(p, Mi)
    double E2 = 0.0;               // E - om
    double P2[3] = {0, 0, 0};      // P - p
    double sig = 0.0;              // sigma_i_cpp(E, P, p, Mi)
    BoostPar bp;                   // parameters of boost(., ., E2, P2)
    double gam = 1.0, beta = 0.0;  // K3E: gam_x, beta_x
    double bh[3] = {0, 0, 0};      // K3E: beta_x_hat
    double bv[3] = {0, 0, 0};      // K3E: beta_x_vec = P2 / E2
    double pms0[3] = {0, 0, 0};    // [v]: om_q1 - om_q2 with q*^2 = kallen(sig, Mx^2, My^2)/(4 sig); v = flavour of the partner
    double ppv1[3] = {0, 0, 0};    // K3B (flavour 1): boost(om, p, E-om, P-p)
    double om_pp1 = 0.0;           // K3B (flavour 1): sqrt(|ppv1|^2 + M1^2)
};

struct K3Const {
    double E, Mtot2, Ecm2, Delta, M1, M2;
    double P[3];
};

inline void k3_momentum(const K3Const& c, const comp E, const std::vector<comp>& total_P, const Group& g, K3Mom& r) {
    using namespace k3_2plus1;
    const double M1 = c.M1, M2 = c.M2;
    const k3_2plus1::Vec3 Pvec = make_vec3(total_P[0], total_P[1], total_P[2]);
    const k3_2plus1::Vec3 pv = make_vec3(g.pc[0], g.pc[1], g.pc[2]);
    r.flav = g.flav;
    r.Mi = flavor_mass(g.flav, M1, M2);
    for (int a = 0; a < 3; ++a) r.p[a] = g.p[a];
    const comp om = omega_vec(pv, r.Mi);
    const comp E2 = E - om;
    const k3_2plus1::Vec3 P2 = Pvec - pv;
    const comp sig = sigma_i_cpp(E, Pvec, pv, r.Mi);
    check_real(om, "omega"); check_real(E2, "E-omega"); check_real(sig, "sigma");
    r.om = om.real(); r.E2 = E2.real(); r.sig = sig.real();
    for (int a = 0; a < 3; ++a) { check_real(P2[a], "P-p"); r.P2[a] = P2[a].real(); }
    {   // parameters of ::boost(., ., E2, P2)
        const comp P2c[3] = {P2[0], P2[1], P2[2]};
        r.bp = make_boost_par(P2c, E2);
    }
    {   // K3E_element: beta_x_vec, gam_x, beta_x, beta_x_hat
        const k3_2plus1::Vec3 bvec = P2 / E2;
        const comp gam = std::sqrt(comp(1.0, 0.0) / (comp(1.0, 0.0) - norm2(bvec)));
        const comp beta = norm3(bvec);
        k3_2plus1::Vec3 bhat = bvec;
        if (std::abs(beta) != 0.0) bhat = bvec / beta;
        check_real(gam, "gamma"); check_real(beta, "beta");
        r.gam = gam.real(); r.beta = beta.real();
        for (int a = 0; a < 3; ++a) { check_real(bvec[a], "beta vec"); check_real(bhat[a], "beta hat"); r.bv[a] = bvec[a].real(); r.bh[a] = bhat[a].real(); }
    }
    // pms_0 = om_q1 - om_q2 with q*^2 = kallentriangle(sig, Mx^2, My^2) / (4 sig); (Mx, My) = (M1, M2) for a kaon partner
    // (v = 1) and (M2, M1) for a pion partner (v = 2); q*^2 + M^2 = (sig + M^2 - M'^2)^2 / (4 sig) >= 0 for a kaon spectator.
    // A PION spectator has no pms_0 here: the oracle computes it there (it is imaginary below the K K threshold, where
    // q*^2 + m_pi^2 < 0) but no branch of K3B_element / K3E_element uses it for i = 2 or j = 2.
    if (g.flav == 1) {
        for (int v = 1; v <= 2; ++v) {
            const double Mx2 = (v == 1 ? M1 * M1 : M2 * M2), My2 = (v == 1 ? M2 * M2 : M1 * M1);
            const comp q2 = kallentriangle(sig, comp(Mx2, 0.0), comp(My2, 0.0)) / (comp(4.0, 0.0) * sig);
            const comp o1 = std::sqrt(q2 + comp(M1 * M1, 0.0)), o2 = std::sqrt(q2 + comp(M2 * M2, 0.0));
            const comp d = o1 - o2;
            check_real(d, "pms_0");
            r.pms0[v] = d.real();
        }
    }
    if (g.flav == 1) {   // K3B_element, i == 1: pp_vec1 = boost(om_pi, pvec, E - om_pi, P2p_vec)
        const k3_2plus1::Vec3 pp = boost_cpp(om, pv, E - om, P2);
        const comp ompp = std::sqrt(norm2(pp) + comp(M1 * M1, 0.0));
        check_real(ompp, "omega_pp1");
        for (int a = 0; a < 3; ++a) { check_real(pp[a], "pp_vec1"); r.ppv1[a] = pp[a].real(); }
        r.om_pp1 = ompp.real();
    }
}

// Block values of the oracle's K3B_element / K3E_element (before the block factor), divided by Mtot2.
// Bv[a][b], Ev[a][b] for a < A.nlm, b < B.nlm; flavour (A,B) in {(1,1),(1,2),(2,2)}.
inline void k3_pair_block(const K3Const& c, const K3Mom& A, const K3Mom& B, const Group& gA, const Group& gB, bool needB, bool needE,
                          double Bv[4][4], double Ev[4][4]) {
    const double E = c.E, M1 = c.M1, M2 = c.M2;
    for (int a = 0; a < gA.nlm; ++a) for (int b = 0; b < gB.nlm; ++b) { Bv[a][b] = 0.0; Ev[a][b] = 0.0; }
    const double omA = A.om, omB = B.om;
    const double* pA = A.p;
    const double* pB = B.p;
    const double* P = c.P;
    const double pk[3] = {pA[0] + pB[0], pA[1] + pB[1], pA[2] + pB[2]};      // pvec + kvec
    const double Pdpk = dot3d(P, pk);                                          // dot3(Pvec, pvec + kvec)
    const double pdk = dot3d(pA, pB);                                          // dot3(pvec, kvec)
    const auto isS = [](int l, int m) { return l == 0 && m == 0; };
    const auto isP = [](int l, int m) { return l == 1 && (m == -1 || m == 0 || m == 1); };
    const int fi = A.flav, fj = B.flav;

    if (fi == 2 && fj == 2) {
        // K3B: val += 2 * (Ecm2 - E*(om_pi + om_kj) + dot3(Pvec, pvec + kvec) + (M2^2 - 4 M1^2)); return val / Mtot2
        // K3E: val += 2 * (M2^2 - om_pi*om_kj + dot3(pvec, kvec)); return val / Mtot2
        for (int a = 0; a < gA.nlm; ++a)
            for (int b = 0; b < gB.nlm; ++b)
                if (isS(gA.ell[a], gA.m[a]) && isS(gB.ell[b], gB.m[b])) {
                    if (needB) Bv[a][b] = (2.0 * (c.Ecm2 - E * (omA + omB) + Pdpk + (M2 * M2 - 4.0 * M1 * M1))) / c.Mtot2;
                    if (needE) Ev[a][b] = (2.0 * (M2 * M2 - omA * omB + pdk)) / c.Mtot2;
                }
        // other (l, m) of pion spectators: no branch matches in the oracle -> 0 (l <= 1 and waves_2 = {0} in production)
        return;
    }
    if (fi == 2 && fj == 1) throw std::logic_error("k3_pair_block: (2,1) is filled by transposition");

    // ---- K3B ----
    if (needB) {
        const double pms0 = A.pms0[1];                      // q*^2 of qst2_i_cpp(E, P, p, M1, M1, M2)
        if (fi == 1 && fj == 1) {
            const double kms0 = B.pms0[1];
            for (int a = 0; a < gA.nlm; ++a)
                for (int b = 0; b < gB.nlm; ++b) {
                    const int la = gA.ell[a], ma = gA.m[a], lb = gB.ell[b], mb = gB.m[b];
                    double val = 0.0;
                    if (isS(la, ma) && isS(lb, mb))
                        // val += E*(om_pi + om_kj) - dot3(Pvec, pvec + kvec) - 6 M1^2 + om_pp1*pms_0 + om_kk1*kms_0
                        val = E * (omA + omB) - Pdpk - 6.0 * M1 * M1 + A.om_pp1 * pms0 + B.om_pp1 * kms0;
                    else if (isP(la, ma) && isS(lb, mb))
                        val = (-2.0 / 3.0) * Yq0(1, ma, A.ppv1);       // -2/3 * y1real(pp_vec1, m_out)
                    else if (isS(la, ma) && isP(lb, mb))
                        val = (-2.0 / 3.0) * Yq0(1, mb, B.ppv1);       // -2/3 * y1real(kk_vec1, m_in)
                    Bv[a][b] = val / c.Mtot2;
                }
        } else {   // (1,2)
            for (int a = 0; a < gA.nlm; ++a)
                for (int b = 0; b < gB.nlm; ++b) {
                    const int la = gA.ell[a], ma = gA.m[a], lb = gB.ell[b], mb = gB.m[b];
                    double val = 0.0;
                    if (isS(la, ma) && isS(lb, mb))
                        // val += E*om_pi - dot3(pvec, Pvec) + Ecm2 - 2*E*om_kj + 2*dot3(Pvec, kvec) + om_pp1*pms_0 + (M2^2 - 7 M1^2)
                        val = E * omA - dot3d(pA, P) + c.Ecm2 - 2.0 * E * omB + 2.0 * dot3d(P, pB) + A.om_pp1 * pms0 + (M2 * M2 - 7.0 * M1 * M1);
                    else if (isP(la, ma) && isS(lb, mb))
                        val = (-2.0 / 3.0) * Yq0(1, ma, A.ppv1);
                    Bv[a][b] = val / c.Mtot2;
                }
        }
    }

    // ---- K3E ----
    if (needE) {
        if (fi == 1 && fj == 2) {
            const double pms0 = A.pms0[2];                  // sig_pi with (Mj^2, Mk^2) = (M2^2, M1^2)
            double ksp[3];
            boost_apply(A.bp, omB, pB, ksp);                // ksp_vec = boost(om_kj, kvec, E2p, P2p_vec)
            const double om_ksp = std::sqrt(dot3d(ksp, ksp) + B.Mi * B.Mi);   // sqrt(norm2(ksp) + Mj^2)
            for (int a = 0; a < gA.nlm; ++a)
                for (int b = 0; b < gB.nlm; ++b) {
                    const int la = gA.ell[a], ma = gA.m[a], lb = gB.ell[b], mb = gB.m[b];
                    double val = 0.0;
                    if (isS(la, ma) && isS(lb, mb))
                        // val += 2 M2^2 - (E - om_pi)*om_kj + dot3(P2p_vec, kvec) + pms_0*om_ksp
                        val = 2.0 * M2 * M2 - A.E2 * omB + dot3d(A.P2, pB) + pms0 * om_ksp;
                    else if (isP(la, ma) && isS(lb, mb))
                        val = (-2.0 / 3.0) * Yq0(1, ma, ksp);
                    Ev[a][b] = val / c.Mtot2;
                }
        } else {   // (1,1)
            const double pms0 = A.pms0[1], kms0 = B.pms0[1];
            double psk[3], ksp[3], P2psk[3], P2ksp[3];
            boost_apply(B.bp, omA, pA, psk);                // psk_vec  = boost(om_pi, pvec, E2k, P2k_vec)
            boost_apply(A.bp, omB, pB, ksp);                // ksp_vec  = boost(om_kj, kvec, E2p, P2p_vec)
            boost_apply(B.bp, A.E2, A.P2, P2psk);           // P2psk_vec = boost(E2p, P2p_vec, E2k, P2k_vec)
            boost_apply(A.bp, B.E2, B.P2, P2ksp);           // P2ksp_vec = boost(E2k, P2k_vec, E2p, P2p_vec)
            const double P2psk_0 = std::sqrt(A.sig + dot3d(P2psk, P2psk));
            const double P2ksp_0 = std::sqrt(B.sig + dot3d(P2ksp, P2ksp));
            const double gam_p = A.gam, gam_k = B.gam, beta_p = A.beta, beta_k = B.beta;
            const double bhat_dot = dot3d(B.bh, A.bh);      // dot3(beta_k_hat, beta_p_hat)
            const double bvec_dot = dot3d(A.bv, B.bv);      // dot3(beta_p_vec, beta_k_vec)
            const double bkhat_bp = dot3d(B.bh, A.bv);      // dot3(beta_k_hat, beta_p_vec)
            const double bphat_bk = dot3d(A.bh, B.bv);      // dot3(beta_p_hat, beta_k_vec)
            double V[3], Vp[3];
            {
                // V  = (-pms_0*gam_p) * (beta_p_vec + dot3(beta_k_hat, beta_p_vec)*(gam_k - 1)*beta_k_hat - gam_k*beta_k_vec)
                // Vp = (-kms_0*gam_k) * (beta_k_vec + dot3(beta_p_hat, beta_k_vec)*(gam_p - 1)*beta_p_hat - gam_p*beta_p_vec)
                const double sV = (bkhat_bp * (gam_k - 1.0)), sVp = (bphat_bk * (gam_p - 1.0));
                const double fV = -pms0 * gam_p, fVp = -kms0 * gam_k;
                for (int a = 0; a < 3; ++a) {
                    V[a] = fV * ((A.bv[a] + sV * B.bh[a]) - gam_k * B.bv[a]);
                    Vp[a] = fVp * ((B.bv[a] + sVp * A.bh[a]) - gam_p * A.bv[a]);
                }
            }
            bool anyPP = false;
            for (int a = 0; a < gA.nlm; ++a) for (int b = 0; b < gB.nlm; ++b) if (isP(gA.ell[a], gA.m[a]) && isP(gB.ell[b], gB.m[b])) anyPP = true;
            double t[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            if (anyPP) {
                // t[Im][Jm]: Im = (I + 2) % 3 maps x, y, z -> local index 2, 0, 1 (m = +1, -1, 0 -> index m + 1)
                for (int I = 0; I < 3; ++I) {
                    const int Im = (I + 2) % 3;
                    t[Im][Im] = -1.0;
                    for (int J = 0; J < 3; ++J) {
                        const int Jm = (J + 2) % 3;
                        t[Im][Jm] += A.bh[I] * B.bh[J] * (gam_k * beta_k * gam_p * beta_p - bhat_dot * (gam_k - 1.0) * (gam_p - 1.0))
                                     - B.bh[I] * B.bh[J] * (gam_k - 1.0)
                                     - A.bh[I] * A.bh[J] * (gam_p - 1.0);
                    }
                }
            }
            double a1[3], a2[3];
            for (int a = 0; a < 3; ++a) { a1[a] = P2ksp[a] + Vp[a]; a2[a] = P2psk[a] + V[a]; }
            for (int a = 0; a < gA.nlm; ++a)
                for (int b = 0; b < gB.nlm; ++b) {
                    const int la = gA.ell[a], ma = gA.m[a], lb = gB.ell[b], mb = gB.m[b];
                    double val = 0.0;
                    if (isS(la, ma) && isS(lb, mb))
                        // val += 2 M2^2 - 0.5*( Ecm2 - E*(om_pi + om_kj) + dot3(Pvec, pvec + kvec) + om_pi*om_kj - dot3(pvec, kvec)
                        //                      - P2psk_0*kms_0 - P2ksp_0*pms_0 + gam_p*gam_k*(1 - bvec_dot)*pms_0*kms_0 )
                        val = 2.0 * M2 * M2
                              - 0.5 * (c.Ecm2 - E * (omA + omB) + Pdpk + omA * omB - pdk
                                       - P2psk_0 * kms0 - P2ksp_0 * pms0 + gam_p * gam_k * (1.0 - bvec_dot) * pms0 * kms0);
                    else if (isP(la, ma) && isS(lb, mb))
                        val = (-1.0 / 3.0) * Yq0(1, ma, a1);          // -1/3 * y1real(P2ksp_vec + Vp, m_out)
                    else if (isS(la, ma) && isP(lb, mb))
                        val = (-1.0 / 3.0) * Yq0(1, mb, a2);          // -1/3 * y1real(P2psk_vec + V, m_in)
                    else if (isP(la, ma) && isP(lb, mb))
                        val = (-2.0 / 3.0) * t[ma + 1][mb + 1];       // -2/3 * t[m_out + 1][m_in + 1]
                    Ev[a][b] = val / c.Mtot2;
                }
        }
    }
}

// The four basis matrices (iso0, iso1, B, E) of K3mat_2plus1 (c = e_i), each N x N real; need[i] selects which to build.
inline void k3_basis_matrices(const qc3_oracle::Kin& kin, const std::vector<Group>& gs, const bool need[4], std::array<MatR, 4>& out) {
    const int N = kin.N;
    const comp Ec = kin.En_c;
    check_real(Ec, "E_lab");
    K3Const c;
    c.E = Ec.real(); c.M1 = kin.par.atmK; c.M2 = kin.par.atmpi;
    {
        const double Mtot_d = 2.0 * c.M1 + c.M2;
        c.Mtot2 = Mtot_d * Mtot_d;
        for (int a = 0; a < 3; ++a) { check_real(kin.total_P[std::size_t(a)], "total_P"); c.P[a] = kin.total_P[std::size_t(a)].real(); }
        const comp Pvec2 = k3_2plus1::norm2(k3_2plus1::make_vec3(kin.total_P[0], kin.total_P[1], kin.total_P[2]));
        const comp Ecm2 = Ec * Ec - Pvec2;
        const comp Delta = Ecm2 / comp(Mtot_d * Mtot_d, 0.0) - comp(1.0, 0.0);   // K3_element_raw
        check_real(Ecm2, "Ecm2"); check_real(Delta, "Delta");
        c.Ecm2 = Ecm2.real(); c.Delta = Delta.real();
    }
    const int ng = int(gs.size());
    std::vector<K3Mom> mom(static_cast<std::size_t>(ng));
    for (int g = 0; g < ng; ++g) k3_momentum(c, Ec, kin.total_P, gs[std::size_t(g)], mom[std::size_t(g)]);

    for (int i = 0; i < 4; ++i) out[std::size_t(i)] = MatR::Zero(N, N);
    const bool nB = need[2], nE = need[3];
    double Bv[4][4], Ev[4][4];
    for (int ia = 0; ia < ng; ++ia) {
        const Group& gA = gs[std::size_t(ia)];
        for (int ib = ia; ib < ng; ++ib) {
            const Group& gB = gs[std::size_t(ib)];
            if (nB || nE) k3_pair_block(c, mom[std::size_t(ia)], mom[std::size_t(ib)], gA, gB, nB, nE, Bv, Ev);
            const double f = k3_2plus1::two_plus_one_block_factor(gA.flav, gB.flav);
            for (int a = 0; a < gA.nlm; ++a)
                for (int b = 0; b < gB.nlm; ++b) {
                    const int r = gA.row0 + a, cc = gB.row0 + b;
                    const bool ss = (gA.ell[a] == 0 && gA.m[a] == 0 && gB.ell[b] == 0 && gB.m[b] == 0);
                    double v[4] = {0.0, 0.0, 0.0, 0.0};
                    if (ss) { v[0] = f * 1.0; v[1] = f * c.Delta; }      // Kiso = {1, 0} -> 1 ; {0, 1} -> Delta
                    if (nB) v[2] = f * Bv[a][b];
                    if (nE) v[3] = f * Ev[a][b];
                    for (int i = 0; i < 4; ++i) {
                        if (!need[i]) continue;
                        out[std::size_t(i)](r, cc) = v[i];
                        if (ia != ib) out[std::size_t(i)](cc, r) = v[i];
                    }
                }
        }
    }
    for (int i = 0; i < 4; ++i)
        if (need[i] && !out[std::size_t(i)].allFinite()) throw std::runtime_error("nonfinite K3 basis matrix (fast pair kernel)");
}

}  // namespace qc3_fast
