// k2inv_forms.hpp : named p-wave forms of K2^{-1} for the oracle (prompts/001 section 2.2).
//
//   coded_v2        : the production function K2inv_EREord2_2plus1_mat, called UNCHANGED (bit-identical).
//   q0_consistent   : new, separately named function below.  Identical to the production code except for the
//                     ell = 1 effective-range bracket, which becomes the Q0-convention form
//                         q^3 cot(delta_1) = -1/a1 + r1 q^2/2 + P1 q^4         (instead of  -1/(a1 q^2) + r1/2 + P1 q^2),
//                     the threshold term q^2 |q| (1-H) is unchanged.  ell = 0 (and ell >= 2) entries are produced by the
//                     unchanged helper K2_inv_ERE_ang_mom.
// Nothing in src/core is modified; the helpers q2psq_star, cutoff_function_1, omega_func are reused.
#pragma once

#include <string>
#include <vector>
#include <stdexcept>

#include "K2_functions_v2.h"

namespace qc3_oracle {

// Same as K2_inv_ERE_ang_mom, except for the ell=1 effective-range bracket.
inline comp K2_inv_ERE_ang_mom_q0_consistent(double eta_i, const std::vector<std::vector<comp>>& scatter_params,
                                             int ell_f, int proj_mf, int ell_i, int proj_mi,
                                             comp sigma_i, double mj, double mk, double epsilon_h) {
    if (ell_f == 1 && ell_i == 1 && proj_mf == proj_mi) {
        const double pi = std::acos(-1.0);
        const comp A = eta_i / (8.0 * pi * std::sqrt(sigma_i));
        const comp C = std::abs(std::sqrt(q2psq_star(sigma_i, mj, mk)));
        const comp D = 1.0 - cutoff_function_1(sigma_i, mj, mk, epsilon_h);
        const comp q2psq = q2psq_star(sigma_i, mj, mk);
        const comp a1 = scatter_params[1][0];
        const comp r1 = scatter_params[1][1];
        const comp P1 = scatter_params[1][2];
        const comp B = -1.0 / a1 + 0.5 * r1 * q2psq + P1 * q2psq * q2psq;   // = q^3 cot(delta_1), Q0 convention
        return A * (B + q2psq * C * D);
    }
    return K2_inv_ERE_ang_mom(eta_i, scatter_params, ell_f, proj_mf, ell_i, proj_mi, sigma_i, mj, mk, epsilon_h);
}

// Diagonal block, same structure as K2inv_EREord2_i_mat (1/(2 omega_k L^3) prefactor), with the q0_consistent p-wave.
inline void K2inv_EREord2_i_mat_q0_consistent(Eigen::MatrixXcd& K2inv, double eta_i,
                                              const std::vector<std::vector<comp>>& scatter_params, comp En,
                                              const std::vector<std::vector<comp>>& klm_config,
                                              const std::vector<comp>& total_P, double mi, double mj, double mk,
                                              double epsilon_h, double L) {
    const int n = int(klm_config[0].size());
    K2inv = Eigen::MatrixXcd::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        const int ell = static_cast<int>(klm_config[3][std::size_t(i)].real());
        const int proj_m = static_cast<int>(klm_config[4][std::size_t(i)].real());
        const comp kx = klm_config[0][std::size_t(i)], ky = klm_config[1][std::size_t(i)], kz = klm_config[2][std::size_t(i)];
        const comp spec_k = std::sqrt(kx * kx + ky * ky + kz * kz);
        const comp Pmk_x = total_P[0] - kx, Pmk_y = total_P[1] - ky, Pmk_z = total_P[2] - kz;
        const comp Pminusk = std::sqrt(Pmk_x * Pmk_x + Pmk_y * Pmk_y + Pmk_z * Pmk_z);
        const comp sig_k = (En - omega_func(spec_k, mi)) * (En - omega_func(spec_k, mi)) - Pminusk * Pminusk;
        const comp val = K2_inv_ERE_ang_mom_q0_consistent(eta_i, scatter_params, ell, proj_m, ell, proj_m, sig_k, mj, mk, epsilon_h);
        const comp omega_k = omega_func(spec_k, mi);
        K2inv(i, i) = (1.0 / (2.0 * omega_k * L * L * L)) * val;
    }
}

// Block diagonal diag(K2inv_1, 2*K2inv_2) as K2inv_EREord2_2plus1_mat.
inline void K2inv_EREord2_2plus1_mat_q0_consistent(Eigen::MatrixXcd& K2inv, double eta_1, double eta_2,
                                                   const std::vector<std::vector<comp>>& scatter_params_1,
                                                   const std::vector<std::vector<comp>>& scatter_params_2, comp En,
                                                   const std::vector<std::vector<comp>>& plm_config,
                                                   const std::vector<std::vector<comp>>& klm_config,
                                                   const std::vector<comp>& total_P, double m1, double m2,
                                                   double epsilon_h, double L) {
    const int s1 = int(plm_config[0].size()), s2 = int(klm_config[0].size());
    Eigen::MatrixXcd K1, K2;
    K2inv_EREord2_i_mat_q0_consistent(K1, eta_1, scatter_params_1, En, plm_config, total_P, m1, m1, m2, epsilon_h, L);
    K2inv_EREord2_i_mat_q0_consistent(K2, eta_2, scatter_params_2, En, klm_config, total_P, m2, m1, m1, epsilon_h, L);
    K2inv = Eigen::MatrixXcd::Zero(s1 + s2, s1 + s2);
    K2inv.topLeftCorner(s1, s1) = K1;
    K2inv.bottomRightCorner(s2, s2) = 2.0 * K2;
}

// Dispatcher used by the oracle.  "coded_v2" is the unchanged production function.
inline void K2inv_2plus1_by_form(const std::string& form, Eigen::MatrixXcd& K2inv, double eta_1, double eta_2,
                                 const std::vector<std::vector<comp>>& sp1, const std::vector<std::vector<comp>>& sp2,
                                 comp En, std::vector<std::vector<comp>>& plm_config,
                                 std::vector<std::vector<comp>>& klm_config, const std::vector<comp>& total_P,
                                 double m1, double m2, double epsilon_h, double L) {
    if (form == "coded_v2")
        K2inv_EREord2_2plus1_mat(K2inv, eta_1, eta_2, sp1, sp2, En, plm_config, klm_config, total_P, m1, m2, epsilon_h, L);
    else if (form == "q0_consistent")
        K2inv_EREord2_2plus1_mat_q0_consistent(K2inv, eta_1, eta_2, sp1, sp2, En, plm_config, klm_config, total_P, m1, m2, epsilon_h, L);
    else
        throw std::invalid_argument("oracle: unknown k2_pwave_form '" + form + "' (coded_v2 | q0_consistent)");
}

}  // namespace qc3_oracle
