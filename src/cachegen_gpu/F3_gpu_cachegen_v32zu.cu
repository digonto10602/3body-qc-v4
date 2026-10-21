// ============================================================================
// F3_gpu_cachegen_v32zu.cu
//
// v32zt GPU cache generator for KKpi F3inv + Vsel caches.
//
// Design requested by Digonto:
//   * Run script takes one Lbyas, or Lbyas_min Lbyas_max [step_size].
//   * For each Lbyas, generate the 5 default irreps:
//       000_A1m, 100_A2, 110_A2, 111_A2, 200_A2
//   * coarseN defaults to 10000.
//   * Build K2inv, F2, G, F3 and Vsel on the GPU.
//   * Invert/solve F3 on the GPU to obtain full F3inv.
//   * Save full F3inv and Vsel, but do NOT build projF3inv and do NOT compute
//     det/eigen/minSV diagnostics in this executable.
//   * Use same-size batched F3inv generation where possible.
//   * Use VRAM-aware chunking.
//   * Print detailed progress and timing so long runs do not look stuck.
//
// Implementation note:
//   This translation unit includes the last working device-resident GPU pipeline
//   (v16 custom batched GETRS) with its main() renamed.  We reuse its GPU kernels,
//   config builders, projection builder, and safe cublas LU/getrs helper.  The new
//   v32zt main below only performs cache generation and metadata writing.
// ============================================================================

#define main v16_device_resident_old_main_disabled
#include "F3_device_resident_true_batch_v16_custom_batched_getrs.cu"
#undef main

#include "projections_gpu_v32y_exact_device.cuh"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "K3_functions_2plus1.hpp"
#include "QC_functions_v2.h"

namespace fvproj_v30j {
struct Convention;
Eigen::MatrixXcd projector_for_convention(
    const std::vector<std::vector<comp>>& plm_config,
    const std::vector<std::vector<int>>& np_config,
    const std::vector<std::vector<comp>>& klm_config,
    const std::vector<std::vector<int>>& nk_config,
    const std::string& irrep,
    const std::vector<comp>& nnP_config,
    int parity,
    const Convention& conv);
}

static constexpr std::uint64_t V32ZT_RECORD_MAGIC = 0x5653325a4f524543ULL; // "VS2ZPREC" mnemonic
#ifdef V32ZU_ALLOBJ_QC_WRITER
static const char* V32ZT_VERSION = "v33d_gpu_component_cachegen_all_ingredients";
#else
static const char* V32ZT_VERSION = "v32zu_gpu_cachegen_F3inv_Vsel_exact_v32y_opt9_concurrent_streams";
#endif

namespace v32zt_cpu_proj {

static int g_vsel_reorthonormalize = 1;
static double g_vsel_orth_tol = 1.0e-10;

static inline std::vector<std::vector<int>> rotations_list()
{
    return {
        { 1,  2,  3}, { 2,  3,  1}, { 3,  1,  2},
        { 1,  3, -2}, { 2, -1,  3}, { 3,  2, -1},
        { 1, -2, -3}, { 2, -3, -1}, { 3, -1, -2},
        { 1, -3,  2}, { 2,  1, -3}, { 3, -2,  1},
        {-1,  2, -3}, {-2,  3, -1}, {-3,  1, -2},
        {-1,  3,  2}, {-2, -1, -3}, {-3,  2,  1},
        {-1, -2,  3}, {-2, -3,  1}, {-3, -1,  2},
        {-1, -3, -2}, {-2,  1,  3}, {-3, -2, -1}
    };
}

static inline std::vector<std::vector<int>> all_perms_of_3(std::vector<int> v)
{
    std::vector<std::vector<int>> out;
    std::sort(v.begin(), v.end());
    do {
        out.push_back(v);
    } while (std::next_permutation(v.begin(), v.end()));
    return out;
}

static inline std::vector<std::vector<int>> Oh_list()
{
    std::vector<std::vector<int>> out;
    auto append_perms = [&](const std::vector<int>& v) {
        const auto perms = all_perms_of_3(v);
        out.insert(out.end(), perms.begin(), perms.end());
    };
    append_perms({ 1,  2,  3});
    append_perms({ 1,  2, -3});
    append_perms({ 1, -2,  3});
    append_perms({-1,  2,  3});
    append_perms({ 1, -2, -3});
    append_perms({-1,  2, -3});
    append_perms({-1, -2,  3});
    append_perms({-1, -2, -3});
    return out;
}

static inline bool is_in_rotations_list(const std::vector<int>& R)
{
    const auto rots = rotations_list();
    return std::find(rots.begin(), rots.end(), R) != rots.end();
}

static inline std::vector<int> cubic_transf(const std::vector<int>& vec, const std::vector<int>& p)
{
    if (vec.size() != 3 || p.size() != 3) throw std::runtime_error("cubic_transf: vec and p must both have size 3");
    std::vector<int> out(3);
    for (int j = 0; j < 3; ++j) {
        const int idx = std::abs(p[j]) - 1;
        const int sgn = (p[j] > 0) ? 1 : -1;
        out[j] = sgn * vec[static_cast<std::size_t>(idx)];
    }
    return out;
}

static inline std::vector<std::vector<int>> little_group(const std::vector<comp>& total_nP)
{
    if (total_nP.size() != 3) throw std::runtime_error("little_group: total_nP must have size 3");
    const int x = static_cast<int>(total_nP[0].real());
    const int y = static_cast<int>(total_nP[1].real());
    const int z = static_cast<int>(total_nP[2].real());
    if (x == 0 && y == 0 && z == 0) return Oh_list();
    if (x == 0 && y == 0) {
        return {{ 1,  2,  3}, {-1,  2,  3}, { 1, -2,  3}, {-1, -2,  3},
                { 2,  1,  3}, {-2,  1,  3}, { 2, -1,  3}, {-2, -1,  3}};
    }
    if (x == y && z == 0) {
        return {{1, 2,  3}, {1, 2, -3}, {2, 1,  3}, {2, 1, -3}};
    }
    if (x == y && y == z) {
        return {{1, 2, 3}, {1, 3, 2}, {2, 1, 3}, {2, 3, 1}, {3, 1, 2}, {3, 2, 1}};
    }
    if (x != 0 && y != 0 && x != y && z == 0) {
        return {{1, 2,  3}, {1, 2, -3}};
    }
    if (x == y && y != z) {
        return {{1, 2, 3}, {2, 1, 3}};
    }
    if (x != 0 && y != 0 && z != 0 && x != y && y != z && z != x) {
        return {{1, 2, 3}};
    }
    throw std::runtime_error("little_group: invalid shell input");
}

static inline std::string conj_class(const std::vector<int>& p)
{
    if (p.size() != 3) throw std::runtime_error("conj_class: p must have size 3");
    const std::vector<int> p_abs = {std::abs(p[0]), std::abs(p[1]), std::abs(p[2])};
    int N_negs = 0;
    int N_correct = 0;
    for (int i = 0; i < 3; ++i) {
        if (p[i] < 0) ++N_negs;
        if (p_abs[static_cast<std::size_t>(i)] == i + 1) ++N_correct;
    }
    if (N_correct == 3) {
        if (N_negs == 0) return "E";
        if (N_negs == 2) return "C4^2";
        if (N_negs == 3) return "i";
        if (N_negs == 1) return "sigma_h";
    } else if (N_correct == 0) {
        return ((N_negs % 2) == 0) ? "C3" : "S6";
    } else if (N_correct == 1) {
        int i_correct = -1;
        for (int i = 0; i < 3; ++i) {
            if (p_abs[static_cast<std::size_t>(i)] == i + 1) { i_correct = i; break; }
        }
        if (i_correct == -1) throw std::runtime_error("conj_class: internal error finding correct index");
        if ((N_negs % 2) == 1) return (p[static_cast<std::size_t>(i_correct)] < 0) ? "C2" : "C4";
        return (p[static_cast<std::size_t>(i_correct)] > 0) ? "sigma_d" : "S4";
    }
    throw std::runtime_error("Error in conj_class: should never reach here");
}

static inline int chi(const std::vector<int>& p, const std::string& I, const std::vector<comp>& Pvec)
{
    if (Pvec.size() != 3) throw std::runtime_error("chi: Pvec must have size 3");
    const std::string cc = conj_class(p);
    const int Px = static_cast<int>(Pvec[0].real());
    const int Py = static_cast<int>(Pvec[1].real());
    const int Pz = static_cast<int>(Pvec[2].real());

    if (Px == 0 && Py == 0 && Pz == 0) {
        if (I == "A1g" || I == "A1") return 1;
        if (I == "A2g" || I == "A2") return (cc == "C2" || cc == "C4" || cc == "sigma_d" || cc == "S4") ? -1 : 1;
        if (I == "Eg" || I == "E") return (cc == "E" || cc == "C4^2" || cc == "i" || cc == "sigma_h") ? 2 : ((cc == "C3" || cc == "S6") ? -1 : 0);
        if (I == "T1g" || I == "T1") return (cc == "E" || cc == "i") ? 3 : ((cc == "C3" || cc == "S6") ? 0 : ((cc == "C4" || cc == "S4") ? 1 : -1));
        if (I == "T2g" || I == "T2") return (cc == "E" || cc == "i") ? 3 : ((cc == "C3" || cc == "S6") ? 0 : ((cc == "C2" || cc == "sigma_d") ? 1 : -1));
        if (I == "A1u") return (cc == "E" || cc == "C3" || cc == "C4^2" || cc == "C4" || cc == "C2") ? 1 : -1;
        if (I == "A2u") return (cc == "E" || cc == "C3" || cc == "C4^2" || cc == "S4" || cc == "sigma_d") ? 1 : -1;
        if (I == "Eu") return (cc == "E" || cc == "C4^2") ? 2 : ((cc == "i" || cc == "sigma_h") ? -2 : ((cc == "S6") ? 1 : ((cc == "C3") ? -1 : 0)));
        if (I == "T1u") return (cc == "E") ? 3 : ((cc == "i") ? -3 : ((cc == "C4" || cc == "sigma_h" || cc == "sigma_d") ? 1 : ((cc == "C2" || cc == "C4^2" || cc == "S4") ? -1 : 0)));
        if (I == "T2u") return (cc == "E") ? 3 : ((cc == "i") ? -3 : ((cc == "C2" || cc == "S4" || cc == "sigma_h") ? 1 : ((cc == "C4" || cc == "C4^2" || cc == "sigma_d") ? -1 : 0)));
        throw std::runtime_error("chi: invalid irrep for Oh");
    }

    if (Px == 0 && Py == 0) {
        if (I == "A1") return 1;
        if (I == "A2") {
            if (cc == "sigma_h" || cc == "sigma_d") return -1;
            if (cc == "E" || cc == "C4" || cc == "C4^2") return 1;
            throw std::runtime_error("chi: " + cc + " not in LG(00a)");
        }
        if (I == "B1") {
            if (cc == "C4" || cc == "sigma_d") return -1;
            if (cc == "E" || cc == "C4^2" || cc == "sigma_h") return 1;
            throw std::runtime_error("chi: " + cc + " not in LG(00a)");
        }
        if (I == "B2") {
            if (cc == "C4" || cc == "sigma_h") return -1;
            if (cc == "E" || cc == "C4^2" || cc == "sigma_d") return 1;
            throw std::runtime_error("chi: " + cc + " not in LG(00a)");
        }
        if (I == "E" || I == "E2") {
            if (cc == "E") return 2;
            if (cc == "C4^2") return -2;
            if (cc == "C4" || cc == "sigma_h" || cc == "sigma_d") return 0;
            throw std::runtime_error("chi: " + cc + " not in LG(00a)");
        }
        throw std::runtime_error("chi: invalid irrep for C4v");
    }

    if (Px == Py && Py != Pz && Pz == 0) {
        if (I == "A1") return 1;
        if (I == "A2") {
            if (cc == "sigma_h" || cc == "sigma_d") return -1;
            if (cc == "E" || cc == "C2") return 1;
            throw std::runtime_error("chi: " + cc + " not in LG(aa0)");
        }
        if (I == "B1") {
            if (cc == "C2" || cc == "sigma_h") return -1;
            if (cc == "E" || cc == "sigma_d") return 1;
            throw std::runtime_error("chi: " + cc + " not in LG(aa0)");
        }
        if (I == "B2") {
            if (cc == "C2" || cc == "sigma_d") return -1;
            if (cc == "E" || cc == "sigma_h") return 1;
            throw std::runtime_error("chi: " + cc + " not in LG(aa0)");
        }
        throw std::runtime_error("chi: invalid irrep for C2v");
    }

    if (Px == Py && Py == Pz) {
        if (I == "A1") return 1;
        if (I == "A2") {
            if (cc == "sigma_d") return -1;
            if (cc == "E" || cc == "C3") return 1;
            throw std::runtime_error("chi: " + cc + " not in LG(aaa)");
        }
        if (I == "E" || I == "E2") {
            if (cc == "E") return 2;
            if (cc == "C3") return -1;
            if (cc == "sigma_d") return 0;
            throw std::runtime_error("chi: " + cc + " not in LG(aaa)");
        }
        throw std::runtime_error("chi: invalid irrep for C3v");
    }

    if ((Px != 0 && Px != Py && Py != 0 && Pz == 0) || (Px == Py && Py != Pz)) {
        if (I == "A1" || I == "A") return 1;
        if (I == "A2" || I == "B") {
            if (cc == "sigma_h" || cc == "sigma_d") return -1;
            if (cc == "E") return 1;
            throw std::runtime_error("chi: " + cc + " not in LG(ab0/aab)");
        }
        throw std::runtime_error("chi: invalid irrep for C2");
    }

    throw std::runtime_error("chi: invalid Pvec input");
}

static inline int cart_axis_from_m_v30m(int m)
{
    if (m == 1) return 0;
    if (m == 0) return 2;
    if (m == -1) return 1;
    return -1;
}

static inline double D_real_element_v30m(int ell, int m_tgt, int m_src, const std::vector<int>& R)
{
    if (ell == 0) return (m_tgt == 0 && m_src == 0) ? 1.0 : 0.0;
    if (ell == 1) {
        const int a_t = static_cast<int>(cart_axis_from_m_v30m(m_tgt));
        const int a_s = static_cast<int>(cart_axis_from_m_v30m(m_src));
        if (a_t < 0 || a_s < 0 || R.size() != 3) return 0.0;
        const int src_axis_for_target = std::abs(R[static_cast<std::size_t>(a_t)]) - 1;
        const int sign = (R[static_cast<std::size_t>(a_t)] > 0) ? 1 : -1;
        return (src_axis_for_target == a_s) ? double(sign) : 0.0;
    }
    return 0.0;
}

static inline Eigen::MatrixXcd S_single_flavor(
    const std::vector<std::vector<comp>>& lm_config,
    const std::vector<std::vector<int>>& n_config,
    const std::vector<int>& R)
{
    const int dim = static_cast<int>(n_config[0].size());
    Eigen::MatrixXcd S = Eigen::MatrixXcd::Zero(dim, dim);
    for (int source = 0; source < dim; ++source) {
        const std::vector<int> n_src = {n_config[0][source], n_config[1][source], n_config[2][source]};
        const std::vector<int> Rn = cubic_transf(n_src, R);
        const int ell_src = static_cast<int>(std::llround(std::real(lm_config[3][source])));
        const int m_src = static_cast<int>(std::llround(std::real(lm_config[4][source])));
        for (int target = 0; target < dim; ++target) {
            const std::vector<int> n_tgt = {n_config[0][target], n_config[1][target], n_config[2][target]};
            if (Rn != n_tgt) continue;
            const int ell_tgt = static_cast<int>(std::llround(std::real(lm_config[3][target])));
            const int m_tgt = static_cast<int>(std::llround(std::real(lm_config[4][target])));
            if (ell_tgt == ell_src && m_tgt == m_src) S(target, source) = comp(1.0, 0.0);
        }
    }
    return S;
}

static inline Eigen::MatrixXcd W_single_flavor(
    const std::vector<std::vector<comp>>& lm_config,
    const std::vector<std::vector<int>>& n_config,
    int parity,
    const std::vector<int>& R)
{
    const int dim = static_cast<int>(n_config[0].size());
    Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(dim, dim);
    for (int source = 0; source < dim; ++source) {
        const std::vector<int> n_src = {n_config[0][source], n_config[1][source], n_config[2][source]};
        const int ell_src = static_cast<int>(std::llround(std::real(lm_config[3][source])));
        const int m_src = static_cast<int>(std::llround(std::real(lm_config[4][source])));
        for (int target = 0; target < dim; ++target) {
            const std::vector<int> n_tgt = {n_config[0][target], n_config[1][target], n_config[2][target]};
            if (n_tgt != n_src) continue;
            const int ell_tgt = static_cast<int>(std::llround(std::real(lm_config[3][target])));
            if (ell_tgt != ell_src) continue;
            const int m_tgt = static_cast<int>(std::llround(std::real(lm_config[4][target])));
            const double pf = (parity == -1 && !is_in_rotations_list(R)) ? -1.0 : 1.0;
            const double d = D_real_element_v30m(ell_src, m_tgt, m_src, R);
            W(target, source) += comp(pf * d, 0.0);
        }
    }
    return W;
}

static inline Eigen::MatrixXcd U_single_flavor(
    const std::vector<std::vector<comp>>& lm_config,
    const std::vector<std::vector<int>>& n_config,
    int parity,
    const std::vector<int>& R)
{
    return S_single_flavor(lm_config, n_config, R) * W_single_flavor(lm_config, n_config, parity, R);
}

static inline Eigen::MatrixXcd U_2plus1(
    const std::vector<std::vector<comp>>& plm_config,
    const std::vector<std::vector<int>>& np_config,
    const std::vector<std::vector<comp>>& klm_config,
    const std::vector<std::vector<int>>& nk_config,
    int parity,
    const std::vector<int>& R)
{
    const Eigen::MatrixXcd U1 = U_single_flavor(plm_config, np_config, parity, R);
    const Eigen::MatrixXcd U2 = U_single_flavor(klm_config, nk_config, parity, R);
    Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(U1.rows() + U2.rows(), U1.cols() + U2.cols());
    U.block(0, 0, U1.rows(), U1.cols()) = U1;
    U.block(U1.rows(), U1.cols(), U2.rows(), U2.cols()) = U2;
    return U;
}

static inline Eigen::MatrixXcd projector_for_convention_v30q_fixed(
    const std::vector<std::vector<comp>>& plm_config,
    const std::vector<std::vector<int>>& np_config,
    const std::vector<std::vector<comp>>& klm_config,
    const std::vector<std::vector<int>>& nk_config,
    const std::string& irrep,
    const std::vector<comp>& nnP_config,
    int parity)
{
    const std::vector<std::vector<int>> LG = little_group(nnP_config);
    const int dI = (irrep == "A1g" || irrep == "A1" || irrep == "A2g" || irrep == "A2") ? 1
                 : ((irrep == "Eg" || irrep == "E" || irrep == "Eu" || irrep == "E2") ? 2
                 : 3);
    const int dim = static_cast<int>(np_config[0].size()) + static_cast<int>(nk_config[0].size());
    Eigen::MatrixXcd P = Eigen::MatrixXcd::Zero(dim, dim);
    for (const auto& R : LG) {
        const double chiR = static_cast<double>(chi(R, irrep, nnP_config));
        P += chiR * U_2plus1(plm_config, np_config, klm_config, nk_config, parity, R);
    }
    if (!LG.empty()) P *= (static_cast<double>(dI) / static_cast<double>(LG.size()));
    return 0.5 * (P + P.adjoint());
}

static inline Eigen::MatrixXcd modified_gram_schmidt_complex(const Eigen::MatrixXcd& V, double tol)
{
    if (V.rows() == 0 || V.cols() == 0) return V;
    Eigen::MatrixXcd Q = Eigen::MatrixXcd::Zero(V.rows(), V.cols());
    for (int j = 0; j < V.cols(); ++j) {
        Eigen::VectorXcd q = V.col(j);
        for (int k = 0; k < j; ++k) {
            q -= Q.col(k) * (Q.col(k).adjoint() * q)(0, 0);
        }
        for (int k = 0; k < j; ++k) {
            q -= Q.col(k) * (Q.col(k).adjoint() * q)(0, 0);
        }
        const double n = q.norm();
        if (!(n > tol) || !std::isfinite(n)) {
            throw std::runtime_error("VSEL_GRAM_SCHMIDT_DEPENDENT_COLUMN");
        }
        Q.col(j) = q / n;
    }
    return Q;
}

static inline void build_projector_from_eigenvectors_near_one(
    const Eigen::MatrixXcd& P_I,
    Eigen::MatrixXcd& Vsel,
    Eigen::MatrixXcd& Pproj,
    double eig_tol = 0.05,
    double norm_tol = 1e-12,
    double proj_tol = 1e-10,
    char debug = 'n')
{
    if (P_I.rows() != P_I.cols()) throw std::runtime_error("build_projector_from_eigenvectors_near_one: P_I must be square.");
    const int N = P_I.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(P_I);
    if (es.info() != Eigen::Success) throw std::runtime_error("build_projector_from_eigenvectors_near_one: eigensolver failed.");
    const Eigen::VectorXd evals = es.eigenvalues();
    const Eigen::MatrixXcd evecs = es.eigenvectors();
    std::vector<int> keep_inds;
    for (int i = 0; i < evals.size(); ++i) {
        if (evals(i) >= 1.0 - eig_tol && evals(i) <= 1.0 + eig_tol) keep_inds.push_back(i);
    }
    const int r = static_cast<int>(keep_inds.size());
    if (r == 0) {
        Vsel = Eigen::MatrixXcd(N, 0);
        Pproj = Eigen::MatrixXcd::Zero(N, N);
        return;
    }
    Vsel.resize(N, r);
    for (int j = 0; j < r; ++j) {
        Vsel.col(j) = evecs.col(keep_inds[static_cast<std::size_t>(j)]);
        const double nrm = Vsel.col(j).norm();
        if (nrm > norm_tol) Vsel.col(j) /= nrm;
    }
    if (g_vsel_reorthonormalize) {
        Vsel = modified_gram_schmidt_complex(Vsel, g_vsel_orth_tol);
    }
    Pproj = Vsel * Vsel.adjoint();
    (void)proj_tol;
    (void)debug;
}

} // namespace v32zt_cpu_proj

namespace fvproj_v30j {
using v32zt_cpu_proj::chi;
using v32zt_cpu_proj::cubic_transf;
using v32zt_cpu_proj::is_in_rotations_list;
using v32zt_cpu_proj::little_group;
}

#include "fv_projector_cartesian_l1_v30m.hpp"  // v4: was ../../k3iso01_swave_cpu_cache_qc_omp_package/source/ (package layout)

#ifndef USE_CONCURRENT_STREAMS
#define USE_CONCURRENT_STREAMS 1
#endif

static constexpr double CONCURRENT_STREAMS_VRAM_SAFETY_FRACTION = 0.85;

// v32zu OPT-1..8: runtime switches keep validated physics/cache paths available while enabling guarded optimizations.
struct V32zuRunOptions
{
    bool use_fused_f2k2 = true;
    bool validate_fused_f2k2 = false;
    bool use_kinematic_precompute = true;
    bool use_stream_pipeline = true;
    bool use_concurrent_streams = true;
    bool use_multi_gpu = false;
    bool use_mixed_precision_f2k2 = false;
    double gpu_memory_safety_fraction = 0.85;
    double concurrent_streams_safety_fraction = CONCURRENT_STREAMS_VRAM_SAFETY_FRACTION;
    int max_batch_energies = 0;
    int device_count = 0;
    std::string irrep_filter_csv;
};

static const char* v32zu_cublas_status_string(cublasStatus_t st)
{
    // v32zu streams: human-readable cuBLAS diagnostics for wave launch failures.
    switch (st) {
        case CUBLAS_STATUS_SUCCESS: return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED: return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE: return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH: return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR: return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED: return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR: return "CUBLAS_STATUS_INTERNAL_ERROR";
        case CUBLAS_STATUS_NOT_SUPPORTED: return "CUBLAS_STATUS_NOT_SUPPORTED";
        case CUBLAS_STATUS_LICENSE_ERROR: return "CUBLAS_STATUS_LICENSE_ERROR";
        default: return "CUBLAS_STATUS_UNKNOWN";
    }
}

#define CHECK_CUDA(call) do { \
    cudaError_t err__ = (call); \
    if (err__ != cudaSuccess) { \
        std::ostringstream os__; \
        os__ << "CUDA error at " << __FILE__ << ":" << __LINE__ \
             << " in " #call " code=" << static_cast<int>(err__) \
             << " (" << cudaGetErrorName(err__) << ") " << cudaGetErrorString(err__); \
        throw std::runtime_error(os__.str()); \
    } \
} while(0)

#define CHECK_CUBLAS(call) do { \
    cublasStatus_t st__ = (call); \
    if (st__ != CUBLAS_STATUS_SUCCESS) { \
        std::ostringstream os__; \
        os__ << "cuBLAS error at " << __FILE__ << ":" << __LINE__ \
             << " in " #call " code=" << static_cast<int>(st__) \
             << " (" << v32zu_cublas_status_string(st__) << ")"; \
        throw std::runtime_error(os__.str()); \
    } \
} while(0)

// v32zu OPT-1: validation metrics are accumulated across chunks without changing cache output.
struct V32zuF2K2ValidationStats
{
    double max_abs_diff_F2 = 0.0;
    double max_rel_diff_F2 = 0.0;
    double max_abs_diff_K2inv = 0.0;
    double max_rel_diff_K2inv = 0.0;
    std::size_t chunks = 0;
};

static V32zuF2K2ValidationStats g_v32zu_f2k2_validation;

static bool v32zu_parse_bool(const std::string& s)
{
    // v32zu config: accept shell-friendly booleans from run scripts.
    return s == "1" || s == "true" || s == "TRUE" || s == "y" || s == "Y" || s == "yes" || s == "YES";
}

struct V32ztStopwatch
{
    std::chrono::high_resolution_clock::time_point t0;
    V32ztStopwatch() : t0(std::chrono::high_resolution_clock::now()) {}
    double elapsed() const
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(t1 - t0).count();
    }
    void reset() { t0 = std::chrono::high_resolution_clock::now(); }
};

struct V32ztTiming
{
    double basis_sec = 0.0;
    double f2k2_sec = 0.0;
    double gf3_build_sec = 0.0;
    double solve_sec = 0.0;
    double vsel_sec = 0.0;
    double write_sec = 0.0;
    double total_sec = 0.0;
    std::size_t rows = 0;
    std::size_t chunks = 0;
    std::size_t bytes_written = 0;
    int max_total_dim = 0;
    int max_vdim = 0;
};

struct V32ztCacheRowPayload
{
    // v32zu streams: host-staged row payload lets concurrent waves serialize binary writes deterministically.
    int grid_i = 0;
    double Ecm = 0.0;
    double En = 0.0;
    int dim1 = 0;
    int dim2 = 0;
    int n = 0;
    int vdim = 0;
    std::vector<cuDoubleComplex> F3inv_colmajor;
    std::vector<cuDoubleComplex> Vsel_colmajor;
};

struct V32ztAllObjRowPayload
{
    // v32zu streams: stage all object matrices for deterministic ordered flush when the all-object build is concurrent.
    int grid_i = 0;
    double Ecm = 0.0;
    Eigen::MatrixXcd plm_config;
    Eigen::MatrixXcd klm_config;
    Eigen::MatrixXcd K2inv;
    Eigen::MatrixXcd F2;
    Eigen::MatrixXcd G;
    Eigen::MatrixXcd F3;
    Eigen::MatrixXcd F3inv;
    Eigen::MatrixXcd Vsel;
    Eigen::MatrixXcd projected_F3inv;
};

struct V32zuChunkTask
{
    // v32zu streams: chunk task is one same-size cuBLAS batch, never mixing matrix dimensions.
    int original_order = 0;
    int group_id = 0;
    int chunk_id = 0;
    int dim1 = 0;
    int dim2 = 0;
    int n = 0;
    std::vector<EnergyPrepassInfo> rows;
    std::size_t estimated_vram_bytes = 0;
};

struct V32zuChunkResult
{
    int original_order = 0;
    int chunk_id = 0;
    int n = 0;
    int dim1 = 0;
    int dim2 = 0;
    V32ztTiming timing;
    std::vector<V32ztCacheRowPayload> rows;
    std::vector<V32ztAllObjRowPayload> allobj_rows;
};

struct GroupStream
{
    // v32zu streams: one stream and one cuBLAS handle exist only while a wave is active.
    int                 group_id = 0;
    int                 n = 0;
    int                 dim1 = 0;
    int                 dim2 = 0;
    int                 batch_size = 0;
    cudaStream_t        stream = nullptr;
    cublasHandle_t      handle = nullptr;
    std::size_t         estimated_vram_bytes = 0;
    bool                launched = false;
    bool                completed = false;
};

struct V32ztIrrepSpec
{
    std::string label;
    std::array<int,3> nnP;
    std::string irrep;
    std::string irrep_tag;
};

static std::vector<V32ztIrrepSpec> v32zt_default_irreps()
{
    return {
        {"000_A1m", {0,0,0}, "A1u", "A1m"},
        {"100_A2",  {0,0,1}, "A2",  "A2"},
        {"110_A2",  {1,1,0}, "A2",  "A2"},
        {"111_A2",  {1,1,1}, "A2",  "A2"},
        {"200_A2",  {0,0,2}, "A2",  "A2"}
    };
}

static std::string v32zt_format_seconds(double s)
{
    if (!std::isfinite(s)) return "nan";
    long long is = static_cast<long long>(std::llround(s));
    long long h = is / 3600;
    long long m = (is % 3600) / 60;
    long long sec = is % 60;
    std::ostringstream os;
    os << h << ":" << std::setw(2) << std::setfill('0') << m
       << ":" << std::setw(2) << std::setfill('0') << sec
       << std::setfill(' ');
    return os.str();
}

static std::string v32zt_ltag(double Lbyas)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(6) << Lbyas;
    std::string s = os.str();
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    for (char& c : s) if (c == '.') c = 'p';
    return s;
}

static std::string v32zt_join_path(const std::string& a, const std::string& b)
{
    if (a.empty()) return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

static void v32zt_write_exact(std::ofstream& f, const void* p, std::size_t n, std::size_t& byte_counter)
{
    f.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n));
    if (!f) throw std::runtime_error("binary cache write failed");
    byte_counter += n;
}

template <class T>
static void v32zt_write_scalar(std::ofstream& f, const T& x, std::size_t& byte_counter)
{
    v32zt_write_exact(f, &x, sizeof(T), byte_counter);
}

static void v32zt_write_cache_record(
    std::ofstream& f,
    int grid_i,
    double Ecm,
    double En,
    int dim1,
    int dim2,
    int n,
    int vdim,
    const std::vector<cuDoubleComplex>& hF3inv_colmajor,
    const std::vector<cuDoubleComplex>& hVsel_colmajor,
    std::size_t& byte_counter)
{
    const std::uint64_t magic = V32ZT_RECORD_MAGIC;
    const std::int32_t i32 = static_cast<std::int32_t>(grid_i);
    const std::int32_t d1 = static_cast<std::int32_t>(dim1);
    const std::int32_t d2 = static_cast<std::int32_t>(dim2);
    const std::int32_t nn = static_cast<std::int32_t>(n);
    const std::int32_t vv = static_cast<std::int32_t>(vdim);

    if (static_cast<int>(hF3inv_colmajor.size()) != n*n)
        throw std::runtime_error("record write: F3inv host size mismatch");
    if (static_cast<int>(hVsel_colmajor.size()) != n*vdim)
        throw std::runtime_error("record write: Vsel host size mismatch");

    v32zt_write_scalar(f, magic, byte_counter);
    v32zt_write_scalar(f, i32, byte_counter);
    v32zt_write_scalar(f, d1, byte_counter);
    v32zt_write_scalar(f, d2, byte_counter);
    v32zt_write_scalar(f, nn, byte_counter);
    v32zt_write_scalar(f, vv, byte_counter);
    v32zt_write_scalar(f, Ecm, byte_counter);
    v32zt_write_scalar(f, En, byte_counter);
    v32zt_write_exact(f, hF3inv_colmajor.data(), hF3inv_colmajor.size()*sizeof(cuDoubleComplex), byte_counter);
    v32zt_write_exact(f, hVsel_colmajor.data(), hVsel_colmajor.size()*sizeof(cuDoubleComplex), byte_counter);
}

template <class T>
static void v32zt_write_scalar(std::ofstream& os, const T& x)
{
    os.write(reinterpret_cast<const char*>(&x), sizeof(T));
    if (!os) throw std::runtime_error("payload write failed");
}

static void v32zt_write_string(std::ofstream& os, const std::string& s)
{
    const std::uint32_t n = static_cast<std::uint32_t>(s.size());
    v32zt_write_scalar(os, n);
    if (n) {
        os.write(s.data(), static_cast<std::streamsize>(n));
        if (!os) throw std::runtime_error("payload write failed");
    }
}

static Eigen::MatrixXcd v32zt_config_to_matrix(const std::vector<std::vector<comp>>& cfg)
{
    const int cols = static_cast<int>(cfg.size());
    int rows = 0;
    for (const auto& c : cfg) rows = std::max(rows, static_cast<int>(c.size()));
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(rows, cols);
    for (int c = 0; c < cols; ++c) {
        for (int r = 0; r < static_cast<int>(cfg[static_cast<std::size_t>(c)].size()); ++r) {
            M(r, c) = cfg[static_cast<std::size_t>(c)][static_cast<std::size_t>(r)];
        }
    }
    return M;
}

static Eigen::MatrixXcd v32zt_complex_buffer_to_matrix(const std::vector<cuDoubleComplex>& buf, int rows, int cols)
{
    Eigen::MatrixXcd M(rows, cols);
    if (static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols) != buf.size())
        throw std::runtime_error("complex buffer size mismatch");
    std::memcpy(M.data(), buf.data(), buf.size() * sizeof(cuDoubleComplex));
    return M;
}

static Eigen::MatrixXcd v32zt_complex_buffer_to_matrix(const cuDoubleComplex* buf, int rows, int cols)
{
    Eigen::MatrixXcd M(rows, cols);
    if (rows < 0 || cols < 0) throw std::runtime_error("complex buffer size mismatch");
    const std::size_t count = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
    if (count) std::memcpy(M.data(), buf, count * sizeof(cuDoubleComplex));
    return M;
}

static std::string v32zt_object_path(
    const std::string& outdir,
    double Lbyas,
    double xi,
    const std::string& irrep,
    int coarseN,
    const std::string& object)
{
    std::ostringstream name;
    name << "v32zu_Lbyas" << v32zt_ltag(Lbyas)
         << "_xi" << v32zt_ltag(xi)
         << "_irrep" << irrep
         << "_coarse" << coarseN
         << "_" << object << "_gpu.bin";
    return outdir + "/" + name.str();
}

static std::string v32zt_sha256_file(const std::string& path)
{
    std::string cmd = "sha256sum '" + path + "'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "";
    char buf[256];
    std::string out;
    if (fgets(buf, sizeof(buf), pipe)) out = buf;
    pclose(pipe);
    const auto sp = out.find(' ');
    return sp == std::string::npos ? out : out.substr(0, sp);
}

static std::string v32zt_now_utc()
{
    const auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream os;
    os << std::put_time(std::gmtime(&t), "%Y-%m-%dT%H:%M:%SZ");
    return os.str();
}

static std::string v32zt_hostname_string()
{
    char buf[256] = {0};
    if (gethostname(buf, sizeof(buf)-1) == 0) return buf;
    return "unknown";
}

static void v32zt_write_matrix_record(
    std::ofstream& os,
    int row,
    double Ecm,
    const Eigen::MatrixXcd& M)
{
    v32zt_write_scalar(os, static_cast<std::int32_t>(row));
    v32zt_write_scalar(os, Ecm);
    v32zt_write_scalar(os, static_cast<std::int32_t>(M.rows()));
    v32zt_write_scalar(os, static_cast<std::int32_t>(M.cols()));
    for (int r = 0; r < M.rows(); ++r) {
        for (int c = 0; c < M.cols(); ++c) {
            const double re = M(r, c).real();
            const double im = M(r, c).imag();
            v32zt_write_scalar(os, re);
            v32zt_write_scalar(os, im);
        }
    }
}

struct V32ztAllObjWriter
{
    struct Entry
    {
        std::string object;
        std::string path;
        std::ofstream os;
    };

    std::vector<Entry> entries;

    Entry& get(const std::string& object)
    {
        for (auto& e : entries) if (e.object == object) return e;
        throw std::runtime_error("missing all-object writer entry: " + object);
    }
};

static void v32zt_write_staged_allobj_rows(
    V32ztAllObjWriter& writer,
    const std::vector<V32ztAllObjRowPayload>& rows,
    V32ztTiming& timing)
{
    // v32zu streams: one ordered flush point keeps all-object outputs deterministic.
    V32ztStopwatch twrite;
    for (const auto& row : rows)
    {
        v32zt_write_matrix_record(writer.get("plm_config").os, row.grid_i, row.Ecm, row.plm_config);
        v32zt_write_matrix_record(writer.get("klm_config").os, row.grid_i, row.Ecm, row.klm_config);
        v32zt_write_matrix_record(writer.get("K2inv").os, row.grid_i, row.Ecm, row.K2inv);
        v32zt_write_matrix_record(writer.get("F2").os, row.grid_i, row.Ecm, row.F2);
        v32zt_write_matrix_record(writer.get("G").os, row.grid_i, row.Ecm, row.G);
        v32zt_write_matrix_record(writer.get("F3").os, row.grid_i, row.Ecm, row.F3);
        v32zt_write_matrix_record(writer.get("F3inv").os, row.grid_i, row.Ecm, row.F3inv);
        v32zt_write_matrix_record(writer.get("Vsel").os, row.grid_i, row.Ecm, row.Vsel);
        v32zt_write_matrix_record(writer.get("projected_F3inv").os, row.grid_i, row.Ecm, row.projected_F3inv);
    }
    for (auto& entry : writer.entries) {
        entry.os.flush();
    }
    timing.write_sec += twrite.elapsed();
}

static V32ztAllObjWriter v32zt_open_allobj_writers(
    const std::string& outdir,
    const V32ztIrrepSpec& spec,
    const PhysicsParams& par,
    int coarseN,
    std::size_t row_count)
{
    const std::vector<std::string> objects = {
        "plm_config", "klm_config", "K2inv", "F2", "G", "F3", "F3inv", "Vsel", "projected_F3inv"
    };
    V32ztAllObjWriter writer;
    writer.entries.reserve(objects.size());
    std::filesystem::create_directories(outdir);
    for (const auto& object : objects) {
        V32ztAllObjWriter::Entry entry;
        entry.object = object;
        entry.path = v32zt_object_path(outdir, par.Lbyas, par.xi, spec.label, coarseN, object);
        entry.os.open(entry.path, std::ios::binary | std::ios::trunc);
        if (!entry.os) throw std::runtime_error("could not open all-object payload: " + entry.path);
        entry.os.write("KKPI_COMPONENT_CACHE_V001", 25);
        v32zt_write_string(entry.os, object);
        v32zt_write_string(entry.os, "gpu_v32zu_component_cachegen");
        v32zt_write_scalar(entry.os, par.Lbyas);
        v32zt_write_scalar(entry.os, par.xi);
        v32zt_write_string(entry.os, spec.label);
        v32zt_write_scalar(entry.os, static_cast<std::int32_t>(coarseN));
        v32zt_write_scalar(entry.os, static_cast<std::int32_t>(row_count));
        writer.entries.push_back(std::move(entry));
    }
    return writer;
}

static void v32zt_close_allobj_writers(V32ztAllObjWriter& writer)
{
    for (auto& entry : writer.entries) {
        if (entry.os.is_open()) entry.os.close();
    }
}

static void v32zt_finalize_allobj_sidecars(
    const V32ztAllObjWriter& writer,
    const V32ztIrrepSpec& spec,
    const PhysicsParams& par,
    double E0,
    double E1,
    int coarseN,
    const V32ztTiming& timing,
    int device_id,
    const cudaDeviceProp& prop,
    bool writer_thread)
{
    for (const auto& entry : writer.entries) {
        const std::string meta_path = entry.path + ".meta.json";
        const std::string sha = v32zt_sha256_file(entry.path);
        const std::uintmax_t bytes = std::filesystem::file_size(entry.path);
        std::ofstream m(meta_path);
        if (!m) throw std::runtime_error("could not open metadata sidecar for writing: " + meta_path);
        m << std::setprecision(17);
        m << "{\n";
        m << "  \"version\": \"" << V32ZT_VERSION << "\",\n";
        m << "  \"parent_version\": \"v33d_cpu_gpu_component_cachegen_package\",\n";
        m << "  \"backend\": \"gpu_v32zu_component_cachegen\",\n";
        m << "  \"source_backend\": \"v32zu_gpu_cachegen_opt1_fusedF2K2_nowarnings_package\",\n";
        m << "  \"component_cache_file_format\": \"KKPI_COMPONENT_CACHE_V001\",\n";
        m << "  \"components_saved_from_gpu_buffers\": true,\n";
        m << "  \"precision_mode\": \"fp64\",\n";
        m << "  \"cache_kind\": \"" << entry.object << "\",\n";
        m << "  \"irrep_label\": \"" << spec.label << "\",\n";
        m << "  \"irrep\": \"" << spec.irrep << "\",\n";
        m << "  \"irrep_tag\": \"" << spec.irrep_tag << "\",\n";
        m << "  \"Lbyas\": " << par.Lbyas << ",\n";
        m << "  \"xi\": " << par.xi << ",\n";
        m << "  \"xival\": " << par.xi << ",\n";
        m << "  \"spatial_L\": " << par.xi * par.Lbyas << ",\n";
        m << "  \"Ecm_min\": " << E0 << ",\n";
        m << "  \"Ecm_max\": " << E1 << ",\n";
        m << "  \"coarseN\": " << coarseN << ",\n";
        m << "  \"waves_vec_1\": [0, 1],\n";
        m << "  \"waves_vec_2\": [0],\n";
        m << "  \"scatter_params\": {\"scatter1_00\": 4.04, \"scatter1_10\": -43.2, \"scatter2_00\": 4.12},\n";
        m << "  \"K3iso0\": " << par.K3iso0 << ",\n";
        m << "  \"K3iso1\": " << par.K3iso1 << ",\n";
        m << "  \"K3B_par\": " << par.K3B_par << ",\n";
        m << "  \"K3E_par\": " << par.K3E_par << ",\n";
        m << "  \"matrix_shape\": \"per-row variable\",\n";
        m << "  \"dtype\": \"complex128\",\n";
        m << "  \"layout\": \"row-major records: energy_index,Ecm,rows,cols,re,im\",\n";
        m << "  \"chunked\": true,\n";
        m << "  \"chunk_index_file\": \"\",\n";
        m << "  \"endianness\": \"little\",\n";
        m << "  \"payload_bytes\": " << bytes << ",\n";
        m << "  \"payload_sha256\": \"" << sha << "\",\n";
        m << "  \"created_utc\": \"" << v32zt_now_utc() << "\",\n";
        m << "  \"hostname\": \"" << v32zt_hostname_string() << "\",\n";
        m << "  \"compiler\": \"nvcc\",\n";
        m << "  \"cuda_version\": \"" << CUDART_VERSION << "\",\n";
        m << "  \"gpu_name\": \"" << prop.name << "\",\n";
        m << "  \"gpu_compute_capability\": \"" << prop.major << "." << prop.minor << "\",\n";
        m << "  \"cuda_streams\": " << 1 << ",\n";
        m << "  \"cuda_batch_E\": " << 0 << ",\n";
        m << "  \"pinned_host_staging\": true,\n";
        m << "  \"async_D2H\": true,\n";
        m << "  \"writer_thread\": " << (writer_thread ? "true" : "false") << ",\n";
        m << "  \"timings_sec\": {\n";
        m << "    \"basis\": " << timing.basis_sec << ",\n";
        m << "    \"F2_K2\": " << timing.f2k2_sec << ",\n";
        m << "    \"G_F3\": " << timing.gf3_build_sec << ",\n";
        m << "    \"solve\": " << timing.solve_sec << ",\n";
        m << "    \"Vsel\": " << timing.vsel_sec << ",\n";
        m << "    \"write\": " << timing.write_sec << ",\n";
        m << "    \"total\": " << timing.total_sec << "\n";
        m << "  },\n";
        m << "  \"memory\": {\"peak_gpu_bytes\": 0, \"peak_cpu_rss_bytes\": 0},\n";
        m << "  \"validation\": {\"validated_against\": \"smoke_only\", \"max_abs_diff\": null, \"max_rel_diff\": null, \"tolerance\": 1e-10}\n";
        m << "}\n";
    }
}

__global__ void v32zt_identity_batch_kernel(cuDoubleComplex* X, int B, int n)
{
    long long tid = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    long long total = static_cast<long long>(B) * n * n;
    if (tid >= total) return;
    long long per = static_cast<long long>(n) * n;
    long long loc = tid % per;
    int row = static_cast<int>(loc % n);
    int col = static_cast<int>(loc / n);
    X[tid] = make_cuDoubleComplex(row == col ? 1.0 : 0.0, 0.0);
}

__global__ void v32zu_fused_F2_K2inv_batch_same_shape_kernel(
    cuDoubleComplex* F2_batch_col, // v32zu OPT-1: output F2, still needed later for F3 construction.
    cuDoubleComplex* K_batch_col,  // v32zu OPT-1: output K2inv, still kept separate so cache physics/layout is unchanged.
    int batchCount,                // v32zu OPT-1: number of same-shape energy rows in this chunk.
    int n,                         // v32zu OPT-1: total matrix dimension for this chunk.
    int size1,                     // v32zu OPT-1: first channel dimension.
    int size2,                     // v32zu OPT-1: second channel dimension.
    const double* En_batch,         // v32zu OPT-1: moving-frame energies for each row.
    const double* p1x, const double* p1y, const double* p1z, // v32zu OPT-1: channel-1 momenta.
    const int* ell1, const int* m1cfg,                       // v32zu OPT-1: channel-1 partial-wave labels.
    const double* p2x, const double* p2y, const double* p2z, // v32zu OPT-1: channel-2 momenta.
    const int* ell2, const int* m2cfg,                       // v32zu OPT-1: channel-2 partial-wave labels.
    double Px, double Py, double Pz,                         // v32zu OPT-1: total momentum components.
    double eta_1, double eta_2,                              // v32zu OPT-1: symmetry factors for K2inv.
    k2gpu::ScatterParamsView sp1,                            // v32zu OPT-1: pair-1 scatter parameters.
    k2gpu::ScatterParamsView sp2,                            // v32zu OPT-1: pair-2 scatter parameters.
    double mK, double mpi,                                   // v32zu OPT-1: masses, same meaning as v32zt.
    double L,                                                // v32zu OPT-1: spatial lattice extent.
    double alpha,                                            // v32zu OPT-1: F2 cutoff alpha.
    double epsilon_h,                                        // v32zu OPT-1: cutoff regulator.
    int max_shell_num,                                       // v32zu OPT-1: F2 zeta-shell cutoff.
    bool Q0norm,                                             // v32zu OPT-1: F2/G normalization flag.
    int erfi_max_terms,                                      // v32zu OPT-1: device erfi max terms.
    double erfi_tol)                                         // v32zu OPT-1: device erfi tolerance.
{
    const long long elems_per_E = static_cast<long long>(n) * n;                         // v32zu OPT-1: elements per batched matrix.
    const int e = static_cast<int>(blockIdx.x);                                           // v32zu OPT-2: one CUDA block owns one fixed-shape energy row.
    if (e >= batchCount) return;                                                          // v32zu OPT-2: guard extra row blocks.
    const double En_real = En_batch[e];                                                   // v32zu OPT-3: stage row energy in a register shared by all element-loop iterations.
    const cuDoubleComplex En_c = make_cuDoubleComplex(En_real, 0.0);                      // v32zu OPT-3: avoid reloading En for every matrix element handled by this thread.
    const int off1 = e * size1;                                                          // v32zu OPT-1: channel-1 basis offset.
    const int off2 = e * size2;                                                          // v32zu OPT-1: channel-2 basis offset.

    for (long long loc_ll = threadIdx.x; loc_ll < elems_per_E; loc_ll += blockDim.x)      // v32zu OPT-2: threads in the row block cooperatively cover n*n elements.
    {
        const int loc = static_cast<int>(loc_ll);                                         // v32zu OPT-2: matrix-local column-major element.
        const int row = loc % n;                                                          // v32zu OPT-2: row index.
        const int col = loc / n;                                                          // v32zu OPT-2: column index.
        const long long tid = static_cast<long long>(e) * elems_per_E + loc_ll;           // v32zu OPT-2: global column-major output location.

        cuDoubleComplex f2_val = make_cuDoubleComplex(0.0, 0.0);                          // v32zu OPT-1: default zero for off-block F2 entries.
        k2gpu::Cx k_val(0.0, 0.0);                                                        // v32zu OPT-1: default zero for off-diagonal K2inv entries.

        if (row < size1 && col < size1)                                                   // v32zu OPT-1: channel-1 block.
        {
            f2_val = f2gpu::F2_ang_mom_d(                                                 // v32zu OPT-1: unchanged F2 physics expression.
                En_c,
                p1x[off1+row], p1y[off1+row], p1z[off1+row], ell1[off1+row], m1cfg[off1+row],
                p1x[off1+col], p1y[off1+col], p1z[off1+col], ell1[off1+col], m1cfg[off1+col],
                Px, Py, Pz,
                mK, mK, mpi,
                L, alpha, epsilon_h,
                max_shell_num, Q0norm, erfi_max_terms, erfi_tol);

            if (row == col)                                                               // v32zu OPT-1: K2inv is diagonal in this basis.
            {
                k2gpu::ConfigView cfg1{size1, p1x+off1, p1y+off1, p1z+off1, ell1+off1, m1cfg+off1}; // v32zu OPT-1: channel-1 K2 view.
                k_val = k2gpu::K2inv_diag_element_gpu(En_real, cfg1, row,                 // v32zu OPT-3: use staged En_real, formula unchanged.
                                                      Px, Py, Pz,
                                                      eta_1, sp1,
                                                      mK, mK, mpi,
                                                      epsilon_h, L);
            }
        }
        else if (row >= size1 && col >= size1)                                           // v32zu OPT-1: channel-2 block.
        {
            const int r2 = row - size1;                                                   // v32zu OPT-1: channel-2 row.
            const int c2 = col - size1;                                                   // v32zu OPT-1: channel-2 column.
            f2_val = f2gpu::F2_ang_mom_d(                                                 // v32zu OPT-1: unchanged F2 physics expression.
                En_c,
                p2x[off2+r2], p2y[off2+r2], p2z[off2+r2], ell2[off2+r2], m2cfg[off2+r2],
                p2x[off2+c2], p2y[off2+c2], p2z[off2+c2], ell2[off2+c2], m2cfg[off2+c2],
                Px, Py, Pz,
                mpi, mK, mK,
                L, alpha, epsilon_h,
                max_shell_num, Q0norm, erfi_max_terms, erfi_tol);

            if (r2 == c2)                                                                 // v32zu OPT-1: K2inv remains diagonal in channel 2.
            {
                k2gpu::ConfigView cfg2{size2, p2x+off2, p2y+off2, p2z+off2, ell2+off2, m2cfg+off2}; // v32zu OPT-1: channel-2 K2 view.
                k_val = 2.0 * k2gpu::K2inv_diag_element_gpu(En_real, cfg2, r2,            // v32zu OPT-3: use staged En_real, formula unchanged.
                                                            Px, Py, Pz,
                                                            eta_2, sp2,
                                                            mpi, mK, mK,
                                                            epsilon_h, L);
            }
        }

        F2_batch_col[tid] = f2_val;                                                       // v32zu OPT-1: write F2 once from fused kernel.
        K_batch_col[tid] = make_cuDoubleComplex(k_val.re, k_val.im);                      // v32zu OPT-1: write K2inv once from fused kernel.
    }
}

static void v32zu_update_diff_stats(
    const std::vector<cuDoubleComplex>& ref,
    const std::vector<cuDoubleComplex>& got,
    double& max_abs,
    double& max_rel)
{
    // v32zu OPT-1 validation: compare column-major complex outputs with FP64 tolerances.
    if (ref.size() != got.size()) throw std::runtime_error("F2/K2 validation size mismatch");
    for (std::size_t i=0; i<ref.size(); ++i) {
        const double dr = cuCreal(ref[i]) - cuCreal(got[i]);
        const double di = cuCimag(ref[i]) - cuCimag(got[i]);
        const double ad = std::hypot(dr, di);
        const double rr = std::hypot(cuCreal(ref[i]), cuCimag(ref[i]));
        max_abs = std::max(max_abs, ad);
        max_rel = std::max(max_rel, ad / std::max(1.0, rr));
    }
}

static std::size_t v32zt_estimate_group_bytes_per_energy(int n)
{
    // v32zu OPT-6-lite: keep the same physics/output but stop using the old 16 MB
    // per-energy pad, which forced tiny batches such as 301 rows for n=47.
    // The active pipeline stores many n*n complex matrices simultaneously:
    // F2,G,K,H,X1,F3,F3inv + work/pointers + Vsel<=n*n.
    const std::size_t z = sizeof(cuDoubleComplex); // v32zu: byte size for one complex matrix element.
    const std::size_t n2 = static_cast<std::size_t>(n) * static_cast<std::size_t>(n); // v32zu: elements in one n x n matrix.
    return 18 * n2 * z + 2 * 1024 * 1024; // v32zu: safer dimension-aware estimate with 2 MB row pad.
}

static std::size_t get_free_vram_bytes()
{
    // v32zu streams: wave packing uses live free VRAM, with allocation-time fallback still authoritative.
    std::size_t free_b = 0, total_b = 0;
    cudaError_t e = cudaMemGetInfo(&free_b, &total_b);
    if (e != cudaSuccess) return 0;
    return free_b;
}

static std::size_t v32zt_process_rss_bytes()
{
    std::ifstream in("/proc/self/status");
    if (!in) return 0;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.rfind("VmRSS:", 0) != 0) continue;
        std::istringstream is(line.substr(6));
        long long kb = 0;
        std::string unit;
        if (is >> kb >> unit)
        {
            if (kb > 0) return static_cast<std::size_t>(kb) * 1024ULL;
        }
    }
    return 0;
}

static std::size_t v32zt_process_swap_bytes()
{
    std::ifstream in("/proc/self/status");
    if (!in) return 0;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.rfind("VmSwap:", 0) != 0) continue;
        std::istringstream is(line.substr(7));
        long long kb = 0;
        std::string unit;
        if (is >> kb >> unit)
        {
            if (kb > 0) return static_cast<std::size_t>(kb) * 1024ULL;
        }
    }
    return 0;
}

static std::size_t v32zt_host_rss_budget_bytes()
{
    const long pages = sysconf(_SC_PHYS_PAGES);
    const long page_sz = sysconf(_SC_PAGE_SIZE);
    if (pages <= 0 || page_sz <= 0) return 0;
    return static_cast<std::size_t>(pages) * static_cast<std::size_t>(page_sz);
}

static std::size_t estimate_group_vram_bytes(const std::vector<EnergyPrepassInfo>& group)
{
    // v32zu streams: include the actual active pipeline buffers plus conservative overhead.
    if (group.empty()) return 0;
    const int B = static_cast<int>(group.size());
    const int dim1 = group[0].dim1;
    const int dim2 = group[0].dim2;
    const int n = group[0].total_dim;
    const std::size_t z = sizeof(cuDoubleComplex);
    const std::size_t d = sizeof(double);
    const std::size_t i = sizeof(int);
    const std::size_t ge = sizeof(ggpu::ConfigEntry);
    const std::size_t ptr = sizeof(cuDoubleComplex*);
    const std::size_t Bsz = static_cast<std::size_t>(B);
    const std::size_t nsz = static_cast<std::size_t>(n);
    const std::size_t dim1sz = static_cast<std::size_t>(dim1);
    const std::size_t dim2sz = static_cast<std::size_t>(dim2);
    const std::size_t n2 = nsz * nsz;
    const std::size_t vdim_bound = nsz;

    std::size_t bytes = 0;
    bytes += Bsz * sizeof(double);                         // dEn / row metadata.
    bytes += Bsz * dim1sz * (3*d + 2*i + ge);              // channel-1 basis labels and G config.
    bytes += Bsz * dim2sz * (3*d + 2*i + ge);              // channel-2 basis labels and G config.
    bytes += Bsz * n2 * z * 8;                             // dF2,dK,dG,dH,dX1,dF3,dF3inv,Vsel-bound staging.
    bytes += Bsz * nsz * vdim_bound * z;                   // dVselBatch worst-case.
    bytes += Bsz * ptr * 8;                                // cuBLAS pointer arrays.
    bytes += Bsz * nsz * sizeof(int) * 2;                  // LU pivots for H and F3 solves.
    bytes += Bsz * sizeof(int) * 2;                        // LU info arrays.
    bytes += 8 * 1024 * 1024;                              // cuBLAS/cuSOLVER/scatter/staging pad.
    bytes = static_cast<std::size_t>(std::ceil(double(bytes) * 1.20)); // alignment and allocator overhead.
    return std::max(bytes, v32zt_estimate_group_bytes_per_energy(n) * Bsz);
}

static double v32zt_baseline_rows_per_second()
{
    // v32zu scheduler calibration: coarseN=1000 / 283.672794974 sec from the latest representative GPU benchmark.
    return 3.525;
}

static double v32zt_predict_task_seconds(const V32zuChunkTask& task)
{
    const double base_rows_per_sec = std::max(1.0e-12, v32zt_baseline_rows_per_second());
    const double baseline_n = 189.0;
    const double size_factor = std::pow(std::max(1.0, double(task.n)) / baseline_n, 2.0);
    return (double(task.rows.size()) / base_rows_per_sec) * size_factor;
}

static bool v32zt_should_use_gpu_lane(
    const V32zuChunkTask& task,
    const V32zuRunOptions& opt,
    std::size_t free_vram_bytes,
    std::size_t rss_bytes)
{
    if (!opt.use_concurrent_streams) return false;
    const std::size_t gpu_budget = static_cast<std::size_t>(
        double(free_vram_bytes) * std::min(0.95, std::max(0.10, opt.concurrent_streams_safety_fraction)));
    const std::size_t host_budget = v32zt_host_rss_budget_bytes();
    const std::size_t host_soft_limit = host_budget > 0 ? static_cast<std::size_t>(double(host_budget) * 0.80) : 0;
    const bool gpu_fit = task.estimated_vram_bytes > 0 && task.estimated_vram_bytes <= gpu_budget;
    const bool host_fit = (host_soft_limit == 0) || (rss_bytes < host_soft_limit);
    const bool worthwhile = task.rows.size() >= 2 && v32zt_predict_task_seconds(task) >= 0.25;
    return gpu_fit && host_fit && worthwhile;
}

static int v32zt_auto_batch_cap(int n, int requested_cap, double safety_fraction)
{
    if (requested_cap > 0) return std::max(1, requested_cap);
    std::size_t free_b = 0, total_b = 0;
    cudaError_t e = cudaMemGetInfo(&free_b, &total_b);
    if (e != cudaSuccess || free_b == 0)
    {
        std::cout << "[v32zt VRAM] cudaMemGetInfo failed; using fallback max_batch=8\n";
        return 8;
    }
    safety_fraction = std::min(0.95, std::max(0.10, safety_fraction));
    double usable = double(free_b) * safety_fraction;
    double per = double(v32zt_estimate_group_bytes_per_energy(n));
    int cap = static_cast<int>(std::floor(usable / std::max(1.0, per)));
    return std::max(1, cap);
}

static void v32zt_write_sidecar_json(
    const std::string& meta_path,
    const V32ztIrrepSpec& spec,
    const PhysicsParams& par,
    double E0,
    double E1,
    int coarseN,
    const std::string& cache_path,
    const V32ztTiming& timing,
    const V32zuRunOptions& opt,
    int device_id,
    bool writer_thread)
{
    std::ofstream m(meta_path);
    if (!m) throw std::runtime_error("could not open metadata sidecar for writing: " + meta_path);
    m << std::setprecision(17);
    m << "{\n";
    m << "  \"version\": \"" << V32ZT_VERSION << "\",\n";
    m << "  \"cache_kind\": \"F3inv_plus_Vsel_gpu_full\",\n";
    m << "  \"cache_file\": \"" << cache_path << "\",\n";
    m << "  \"projection_source\": \"k3iso01_swave_cpu_cache_qc_omp_package/source/generate_swave_cpu_cache_qc.cpp canonical CPU Vsel\",\n";
    m << "  \"basis_ordering\": \"config_maker_4_momentum_first\",\n";
    m << "  \"irrep_label\": \"" << spec.label << "\",\n";
    m << "  \"irrep\": \"" << spec.irrep << "\",\n";
    m << "  \"irrep_tag\": \"" << spec.irrep_tag << "\",\n";
    m << "  \"nnP\": [" << spec.nnP[0] << ", " << spec.nnP[1] << ", " << spec.nnP[2] << "],\n";
    m << "  \"Lbyas\": " << par.Lbyas << ",\n";
    m << "  \"xi\": " << par.xi << ",\n";
    m << "  \"spatial_L\": " << par.xi * par.Lbyas << ",\n";
    m << "  \"Ecm_min\": " << E0 << ",\n";
    m << "  \"Ecm_max\": " << E1 << ",\n";
    m << "  \"coarseN\": " << coarseN << ",\n";
    m << "  \"waves_vec_1\": [0, 1],\n";
    m << "  \"waves_vec_2\": [0],\n";
    m << "  \"scatter1_00\": 4.04,\n";
    m << "  \"scatter1_10\": -43.2,\n";
    m << "  \"scatter2_00\": 4.12,\n";
    m << "  \"masses\": {\"atmK\": " << par.atmK << ", \"atmpi\": " << par.atmpi << "},\n";
    m << "  \"alpha\": " << par.alpha << ",\n";
    m << "  \"epsilon_h\": " << par.epsilon_h << ",\n";
    m << "  \"max_shell_num\": " << par.max_shell_num << ",\n";
    m << "  \"Q0norm\": " << (par.Q0norm ? "true" : "false") << ",\n";
    m << "  \"device_id\": " << device_id << ",\n";
    m << "  \"optimization_flags\": {\n";                                                // v32zu OPT diagnostics: additive metadata only, reader layout unchanged.
    m << "    \"use_fused_f2k2\": " << (opt.use_fused_f2k2 ? "true" : "false") << ",\n";
    m << "    \"validate_fused_f2k2\": " << (opt.validate_fused_f2k2 ? "true" : "false") << ",\n";
    m << "    \"use_kinematic_precompute\": " << (opt.use_kinematic_precompute ? "true" : "false") << ",\n";
    m << "    \"use_stream_pipeline_requested\": " << (opt.use_stream_pipeline ? "true" : "false") << ",\n";
    m << "    \"use_stream_pipeline_enabled\": " << (opt.use_concurrent_streams ? "true" : "false") << ",\n";
    m << "    \"use_concurrent_streams\": " << (opt.use_concurrent_streams ? "true" : "false") << ",\n";
    m << "    \"concurrent_streams_safety_fraction\": " << opt.concurrent_streams_safety_fraction << ",\n";
    m << "    \"use_multi_gpu\": " << (opt.use_multi_gpu ? "true" : "false") << ",\n";
    m << "    \"use_mixed_precision_f2k2\": " << (opt.use_mixed_precision_f2k2 ? "true" : "false") << "\n";
    m << "  },\n";
    m << "  \"writer_thread\": " << (writer_thread ? "true" : "false") << ",\n";
    m << "  \"rows_written\": " << timing.rows << ",\n";
    m << "  \"chunks\": " << timing.chunks << ",\n";
    m << "  \"max_total_dim\": " << timing.max_total_dim << ",\n";
    m << "  \"max_vdim\": " << timing.max_vdim << ",\n";
    m << "  \"bytes_written_binary_payload\": " << timing.bytes_written << ",\n";
    m << "  \"timing_seconds\": {\n";
    m << "    \"basis_construction_time\": " << timing.basis_sec << ",\n";
    m << "    \"F2_K2_time\": " << timing.f2k2_sec << ",\n";
    m << "    \"G_F3_matrix_build_time\": " << timing.gf3_build_sec << ",\n";
    m << "    \"matrix_inverse_solve_time\": " << timing.solve_sec << ",\n";
    m << "    \"Vsel_construction_time\": " << timing.vsel_sec << ",\n";
    m << "    \"binary_write_time\": " << timing.write_sec << ",\n";
    m << "    \"total_irrep_time\": " << timing.total_sec << "\n";
    m << "  }\n";
    m << "}\n";
}


static EnergyPrepassInfo v32zt_energy_shape_prepass_momentum_first(
    int idx,
    double Ecm,
    const std::array<int,3>& nnP,
    const PhysicsParams& par,
    char debug)
{
    EnergyPrepassInfo info;
    info.i = idx;
    info.Ecm = Ecm;

    const double L = par.xi * par.Lbyas;
    const double pi = std::acos(-1.0);
    const double twopibyL = 2.0*pi/L;

    std::vector<comp> total_P(3);
    total_P[0] = comp(twopibyL * nnP[0], 0.0);
    total_P[1] = comp(twopibyL * nnP[1], 0.0);
    total_P[2] = comp(twopibyL * nnP[2], 0.0);

    const double P2 = std::norm(total_P[0]) + std::norm(total_P[1]) + std::norm(total_P[2]);
    const double En_real = std::sqrt(Ecm*Ecm + P2);
    info.En = En_real;

    std::vector<int> waves_vec_1 = {0,1};
    std::vector<int> waves_vec_2 = {0};

    std::vector<std::vector<comp>> plm_config(5), klm_config(5);
    std::vector<std::vector<int>> np_config(5), nk_config(5);

    // v32w/v32x-corrected basis ordering: momentum first, then ell/m.
    config_maker_4_momentum_first(plm_config, np_config, waves_vec_1,
                                  comp(En_real,0.0), total_P,
                                  par.atmK, par.atmK, par.atmpi, L,
                                  par.epsilon_h, par.max_shell_num, par.tolerance);
    config_maker_4_momentum_first(klm_config, nk_config, waves_vec_2,
                                  comp(En_real,0.0), total_P,
                                  par.atmpi, par.atmK, par.atmK, L,
                                  par.epsilon_h, par.max_shell_num, par.tolerance);

    info.dim1 = static_cast<int>(plm_config[0].size());
    info.dim2 = static_cast<int>(klm_config[0].size());
    info.total_dim = info.dim1 + info.dim2;

    if (debug == 'y') {
        std::cout << "[v32zt prepass momentum-first] i=" << idx
                  << " Ecm=" << std::setprecision(17) << Ecm
                  << " dim1=" << info.dim1
                  << " dim2=" << info.dim2
                  << " total_dim=" << info.total_dim << "\n";
    }
    return info;
}

static void v32zu_cublas_lu_solve_inplace_batched_same_n_stream(
    cublasHandle_t handle,
    cudaStream_t stream,
    int n,
    int nrhs,
    std::vector<cuDoubleComplex*>& A_lu_ptrs,
    std::vector<cuDoubleComplex*>& B_inout_ptrs,
    const char* label,
    char debug)
{
    // v32zu streams: same-n batched LU on the caller's stream; custom GETRS avoids prior cuBLAS getrs crash.
    const int batch_count = static_cast<int>(A_lu_ptrs.size());
    if (batch_count != static_cast<int>(B_inout_ptrs.size()))
        throw std::runtime_error("streamed batched LU: A/B batch size mismatch");
    if (batch_count <= 0) return;
    if (batch_count == 1)
    {
        // v32zu streams: preserve the validated single-matrix cuSOLVER path for correctness.
        static std::mutex single_lu_mutex;
        CHECK_CUDA(cudaStreamSynchronize(stream));
        std::lock_guard<std::mutex> lock(single_lu_mutex);
        device_lu_solve_inplace(handle, n, nrhs, A_lu_ptrs[0], B_inout_ptrs[0]);
        return;
    }

    DevicePtr<cuDoubleComplex*> dA(static_cast<size_t>(batch_count));
    DevicePtr<cuDoubleComplex*> dB(static_cast<size_t>(batch_count));
    DevicePtr<int> piv(static_cast<size_t>(n) * static_cast<size_t>(batch_count));
    DevicePtr<int> info(static_cast<size_t>(batch_count));

    CHECK_CUDA(cudaMemcpyAsync(dA.p, A_lu_ptrs.data(), sizeof(cuDoubleComplex*)*static_cast<size_t>(batch_count), cudaMemcpyHostToDevice, stream));
    CHECK_CUDA(cudaMemcpyAsync(dB.p, B_inout_ptrs.data(), sizeof(cuDoubleComplex*)*static_cast<size_t>(batch_count), cudaMemcpyHostToDevice, stream));

    if (debug == 'y')
    {
        std::cout << "[v32zu streamed LU] label=" << label
                  << " n=" << n << " nrhs=" << nrhs
                  << " batch=" << batch_count << "\n";
    }

    CHECK_CUBLAS(cublasZgetrfBatched(handle, n, dA.p, n, piv.p, info.p, batch_count));

    std::vector<int> hinfo(static_cast<size_t>(batch_count), 0);
    CHECK_CUDA(cudaMemcpyAsync(hinfo.data(), info.p, sizeof(int)*hinfo.size(), cudaMemcpyDeviceToHost, stream));
    CHECK_CUDA(cudaStreamSynchronize(stream));
    int singular = 0;
    for (int b=0; b<batch_count; ++b)
    {
        const int inf = hinfo[static_cast<size_t>(b)];
        if (inf != 0)
        {
            ++singular;
            std::cout << "[cuBLAS singular info] label=" << label
                      << " batch_local=" << b
                      << " n=" << n
                      << " info=" << inf << "\n";
        }
    }
    if (singular > 0)
    {
        std::cout << "[cuBLAS singular info] label=" << label
                  << " singular_rows=" << singular
                  << " batch=" << batch_count
                  << " note=continuing per requested diagnostics\n";
    }

    dim3 grid(static_cast<unsigned int>(batch_count), static_cast<unsigned int>(nrhs), 1);
    v16_getrs_from_getrf_batched_kernel<<<grid, 1, 0, stream>>>(
        const_cast<const cuDoubleComplex* const*>(dA.p),
        dB.p,
        piv.p,
        n,
        nrhs,
        batch_count);
    CHECK_CUDA(cudaGetLastError());
    CHECK_CUDA(cudaStreamSynchronize(stream));
}

static void v32zt_process_same_shape_chunk_cache(
    const std::vector<EnergyPrepassInfo>& group,
    const V32ztIrrepSpec& spec,
    const PhysicsParams& par,
    std::ofstream& fout,
    V32ztTiming& timing,
    int chunk_id,
    int chunk_total,
    std::size_t global_done_before,
    std::size_t global_total,
    const V32zuRunOptions& opt,
    char debug,
    std::vector<V32ztCacheRowPayload>* staged_rows = nullptr,
    std::vector<V32ztAllObjRowPayload>* staged_allobj_rows = nullptr,
    V32ztAllObjWriter* allobj = nullptr,
    cudaStream_t stream = nullptr,
    cublasHandle_t external_cublas = nullptr,
    std::shared_future<void> launch_gate = std::shared_future<void>(),
    std::promise<void>* ready_promise = nullptr)
{
    if (group.empty()) return;
    const int B = static_cast<int>(group.size());
    const int dim1 = group[0].dim1;
    const int dim2 = group[0].dim2;
    const int n = group[0].total_dim;
    if (n <= 0) return;
    const bool streamed_mode = (stream != nullptr);
    auto sync_point = [&]() {
        // v32zu streams: concurrent chunks synchronize only their own stream.
        if (streamed_mode) CHECK_CUDA(cudaStreamSynchronize(stream));
        else V11_CUDA_CHECK(cudaDeviceSynchronize());
    };

    const double pi = std::acos(-1.0);
    const double L = par.xi * par.Lbyas;
    const double twopibyL = 2.0*pi/L;
    std::vector<comp> total_P(3);
    total_P[0] = comp(twopibyL * spec.nnP[0], 0.0);
    total_P[1] = comp(twopibyL * spec.nnP[1], 0.0);
    total_P[2] = comp(twopibyL * spec.nnP[2], 0.0);

    std::cout << "[v32zu chunk] Lbyas=" << par.Lbyas
              << " irrep=" << spec.label
              << " chunk=" << chunk_id << "/" << chunk_total
              << " batch_rows=" << B
              << " dim1=" << dim1 << " dim2=" << dim2 << " n=" << n
              << " global_done=" << global_done_before << "/" << global_total << "\n";

    V32ztStopwatch t_total;
    V32ztStopwatch t;

    // ---------------- CPU/OpenMP config rebuild for contiguous batch upload ----------------
    std::vector<std::vector<std::vector<comp>>> plm_all(static_cast<std::size_t>(B));
    std::vector<std::vector<std::vector<comp>>> klm_all(static_cast<std::size_t>(B));
    std::vector<std::vector<std::vector<int>>> np_all(static_cast<std::size_t>(B));
    std::vector<std::vector<std::vector<int>>> nk_all(static_cast<std::size_t>(B));
    std::vector<int> waves_vec_1 = {0,1};
    std::vector<int> waves_vec_2 = {0};

    #pragma omp parallel for schedule(dynamic)
    for (int e=0; e<B; ++e) {
        plm_all[static_cast<std::size_t>(e)] = std::vector<std::vector<comp>>(5);
        klm_all[static_cast<std::size_t>(e)] = std::vector<std::vector<comp>>(5);
        np_all[static_cast<std::size_t>(e)]  = std::vector<std::vector<int>>(5);
        nk_all[static_cast<std::size_t>(e)]  = std::vector<std::vector<int>>(5);
        const double En_real = group[static_cast<std::size_t>(e)].En;
        config_maker_4_momentum_first(plm_all[static_cast<std::size_t>(e)], np_all[static_cast<std::size_t>(e)], waves_vec_1,
                       comp(En_real,0.0), total_P,
                       par.atmK, par.atmK, par.atmpi, L, par.epsilon_h,
                       par.max_shell_num, par.tolerance);
        config_maker_4_momentum_first(klm_all[static_cast<std::size_t>(e)], nk_all[static_cast<std::size_t>(e)], waves_vec_2,
                       comp(En_real,0.0), total_P,
                       par.atmpi, par.atmK, par.atmK, L, par.epsilon_h,
                       par.max_shell_num, par.tolerance);
    }
    for (int e=0; e<B; ++e) {
        if ((int)plm_all[e][0].size() != dim1 || (int)klm_all[e][0].size() != dim2) {
            throw std::runtime_error("v32zt chunk rejected: rebuilt config dimensions differ from prepass");
        }
    }

    std::vector<double> hEn(static_cast<std::size_t>(B));
    std::vector<double> p1x(static_cast<std::size_t>(B)*dim1), p1y(p1x.size()), p1z(p1x.size());
    std::vector<int>    e1l(static_cast<std::size_t>(B)*dim1), e1m(e1l.size());
    std::vector<double> p2x(static_cast<std::size_t>(B)*dim2), p2y(p2x.size()), p2z(p2x.size());
    std::vector<int>    e2l(static_cast<std::size_t>(B)*dim2), e2m(e2l.size());
    std::vector<ggpu::ConfigEntry> g1(static_cast<std::size_t>(B)*dim1), g2(static_cast<std::size_t>(B)*dim2);

    for (int e=0; e<B; ++e) {
        hEn[static_cast<std::size_t>(e)] = group[static_cast<std::size_t>(e)].En;
        for (int i=0; i<dim1; ++i) {
            const std::size_t off = static_cast<std::size_t>(e)*dim1 + i;
            p1x[off]=plm_all[e][0][i].real(); p1y[off]=plm_all[e][1][i].real(); p1z[off]=plm_all[e][2][i].real();
            e1l[off]=(int)std::llround(plm_all[e][3][i].real()); e1m[off]=(int)std::llround(plm_all[e][4][i].real());
            g1[off] = ggpu::ConfigEntry{p1x[off], p1y[off], p1z[off], e1l[off], e1m[off]};
        }
        for (int i=0; i<dim2; ++i) {
            const std::size_t off = static_cast<std::size_t>(e)*dim2 + i;
            p2x[off]=klm_all[e][0][i].real(); p2y[off]=klm_all[e][1][i].real(); p2z[off]=klm_all[e][2][i].real();
            e2l[off]=(int)std::llround(klm_all[e][3][i].real()); e2m[off]=(int)std::llround(klm_all[e][4][i].real());
            g2[off] = ggpu::ConfigEntry{p2x[off], p2y[off], p2z[off], e2l[off], e2m[off]};
        }
    }

    DevicePtr<double> dEn(static_cast<std::size_t>(B));
    DevicePtr<double> dp1x(p1x.size()), dp1y(p1y.size()), dp1z(p1z.size());
    DevicePtr<int>    de1l(e1l.size()), de1m(e1m.size());
    DevicePtr<double> dp2x(p2x.size()), dp2y(p2y.size()), dp2z(p2z.size());
    DevicePtr<int>    de2l(e2l.size()), de2m(e2m.size());
    DevicePtr<ggpu::ConfigEntry> dg1(g1.size()), dg2(g2.size());

    V11_CUDA_CHECK(cudaMemcpy(dEn.p, hEn.data(), sizeof(double)*hEn.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(dp1x.p, p1x.data(), sizeof(double)*p1x.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(dp1y.p, p1y.data(), sizeof(double)*p1y.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(dp1z.p, p1z.data(), sizeof(double)*p1z.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(de1l.p, e1l.data(), sizeof(int)*e1l.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(de1m.p, e1m.data(), sizeof(int)*e1m.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(dp2x.p, p2x.data(), sizeof(double)*p2x.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(dp2y.p, p2y.data(), sizeof(double)*p2y.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(dp2z.p, p2z.data(), sizeof(double)*p2z.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(de2l.p, e2l.data(), sizeof(int)*e2l.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(de2m.p, e2m.data(), sizeof(int)*e2m.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(dg1.p, g1.data(), sizeof(ggpu::ConfigEntry)*g1.size(), cudaMemcpyHostToDevice));
    V11_CUDA_CHECK(cudaMemcpy(dg2.p, g2.data(), sizeof(ggpu::ConfigEntry)*g2.size(), cudaMemcpyHostToDevice));

    double basis_sec = t.elapsed();
    timing.basis_sec += basis_sec;
    std::cout << "[timing chunk] basis construction time = " << basis_sec << " sec\n";

    // ---------------- Vsel on GPU ----------------
    t.reset();
    std::vector<int> vdims(static_cast<std::size_t>(B), 0);
    std::vector<cuDoubleComplex*> tempV(static_cast<std::size_t>(B), nullptr);
    int common_vdim = -1;
    for (int e=0; e<B; ++e) {
        v32zt_vsel_gpu::build_projection_Vsel_device_exact_v32y(tempV[static_cast<std::size_t>(e)], vdims[static_cast<std::size_t>(e)], n, dim1,
                                     np_all[static_cast<std::size_t>(e)], nk_all[static_cast<std::size_t>(e)], spec.irrep,
                                     spec.nnP[0], spec.nnP[1], spec.nnP[2], par.parity, 0.05, 1.0e-14, 256, debug);
        if (common_vdim < 0) common_vdim = vdims[static_cast<std::size_t>(e)];
        if (vdims[static_cast<std::size_t>(e)] != common_vdim) {
            for (auto p : tempV) if (p) cudaFree(p);
            throw std::runtime_error("v32zt chunk rejected: vdim is not uniform inside same-shape chunk");
        }
    }
    if (common_vdim <= 0) throw std::runtime_error("v32zt chunk rejected: common_vdim <= 0");
    const int vdim = common_vdim;
    DevicePtr<cuDoubleComplex> dVselBatch(static_cast<std::size_t>(B)*n*vdim);
    for (int e=0; e<B; ++e) {
        V11_CUDA_CHECK(cudaMemcpy(dVselBatch.p + static_cast<std::size_t>(e)*n*vdim,
                                  tempV[static_cast<std::size_t>(e)],
                                  sizeof(cuDoubleComplex)*static_cast<std::size_t>(n)*vdim,
                                  cudaMemcpyDeviceToDevice));
        cudaFree(tempV[static_cast<std::size_t>(e)]);
    }
    V11_CUDA_CHECK(cudaDeviceSynchronize());
    double vsel_sec = t.elapsed();
    timing.vsel_sec += vsel_sec;
    timing.max_vdim = std::max(timing.max_vdim, vdim);
    std::cout << "[timing chunk] Vsel construction time = " << vsel_sec << " sec, vdim=" << vdim << "\n";
    if (streamed_mode && ready_promise)
    {
        // v32zu streams: all wave tasks stop here so Vsel's device-wide sync cannot drain another group's streamed kernels.
        ready_promise->set_value();
        launch_gate.wait();
    }

    // ---------------- Build F2/G/K2/H batches ----------------
    const long long elems_per_E = static_cast<long long>(n)*n;
    const long long total_elems = elems_per_E * B;
    const int threads = 256;                                                             // v32zu OPT-2: row-blocked fused kernel uses 256 threads per energy row.
    const int blocks = static_cast<int>((total_elems + threads - 1)/threads);            // v32zu: flat launch geometry remains for legacy v14 kernels.
    const int fused_blocks = B;                                                          // v32zu OPT-2: one fused F2/K2 block per same-shape energy row.

    DevicePtr<cuDoubleComplex> dF2(static_cast<std::size_t>(total_elems));
    DevicePtr<cuDoubleComplex> dG(static_cast<std::size_t>(total_elems));
    DevicePtr<cuDoubleComplex> dK(static_cast<std::size_t>(total_elems));
    DevicePtr<cuDoubleComplex> dH(static_cast<std::size_t>(total_elems));

    t.reset();
    k2gpu::DeviceScatterParams dSP1; dSP1.upload(make_scatter_params_1()); // v32zu OPT-1: upload scatter params before fused F2/K2 launch.
    k2gpu::DeviceScatterParams dSP2; dSP2.upload(make_scatter_params_2()); // v32zu OPT-1: upload scatter params before fused F2/K2 launch.
    if (opt.use_fused_f2k2) {                                                 // v32zu OPT-1: fused path remains default but is runtime-disablable.
        v32zu_fused_F2_K2inv_batch_same_shape_kernel<<<fused_blocks, threads, 0, stream>>>(// v32zu OPT-9: launch on group stream when concurrent.
            dF2.p, dK.p, B, n, dim1, dim2, dEn.p,                              // v32zu OPT-1: pass both outputs so downstream F3 code is byte-layout unchanged.
            dp1x.p, dp1y.p, dp1z.p, de1l.p, de1m.p,                             // v32zu OPT-1: channel-1 basis arrays.
            dp2x.p, dp2y.p, dp2z.p, de2l.p, de2m.p,                             // v32zu OPT-1: channel-2 basis arrays.
            total_P[0].real(), total_P[1].real(), total_P[2].real(),             // v32zu OPT-1: total momentum.
            par.eta_1, par.eta_2, dSP1.view(), dSP2.view(),                      // v32zu OPT-1: K2inv parameters.
            par.atmK, par.atmpi, L, par.alpha, par.epsilon_h,                    // v32zu OPT-1: masses and F2 parameters.
            (int)par.max_shell_num, par.Q0norm, 100, 1.0e-15);                   // v32zu OPT-1: unchanged F2 summation controls.
        V11_CUDA_CHECK(cudaGetLastError());                                      // v32zu OPT-1: check fused kernel launch.
    } else {
        v14_build_F2_batch_same_shape_kernel<<<blocks, threads, 0, stream>>>(     // v32zu OPT-9: launch on group stream when concurrent.
            dF2.p, B, n, dim1, dim2, dEn.p,
            dp1x.p, dp1y.p, dp1z.p, de1l.p, de1m.p,
            dp2x.p, dp2y.p, dp2z.p, de2l.p, de2m.p,
            total_P[0].real(), total_P[1].real(), total_P[2].real(),
            par.atmK, par.atmpi, L, par.alpha, par.epsilon_h,
            (int)par.max_shell_num, par.Q0norm, 100, 1.0e-15);
        V11_CUDA_CHECK(cudaGetLastError());
        v14_build_K2_batch_same_shape_kernel<<<blocks, threads, 0, stream>>>(     // v32zu OPT-9: launch on group stream when concurrent.
            dK.p, B, n, dim1, dim2, dEn.p,
            dp1x.p, dp1y.p, dp1z.p, de1l.p, de1m.p,
            dp2x.p, dp2y.p, dp2z.p, de2l.p, de2m.p,
            total_P[0].real(), total_P[1].real(), total_P[2].real(),
            par.eta_1, par.eta_2, dSP1.view(), dSP2.view(),
            par.atmK, par.atmpi, par.epsilon_h, L);
        V11_CUDA_CHECK(cudaGetLastError());
    }
    sync_point();                                                            // v32zu OPT-9: stream-local timing in concurrent mode.
    double f2k2_sec = t.elapsed();
    timing.f2k2_sec += f2k2_sec;
    std::cout << "[F2K2] mode=" << (opt.use_fused_f2k2 ? "fused_row_blocked" : "separate_v14")
              << " blocks=" << (opt.use_fused_f2k2 ? fused_blocks : blocks)
              << " threads=" << threads << "\n";
    std::cout << "[timing chunk] fused F2/K2 time = " << f2k2_sec << " sec\n";

    if (opt.validate_fused_f2k2) {                                             // v32zu OPT-1 validation: run old and new kernels on the same small chunk.
        DevicePtr<cuDoubleComplex> dF2_ref(static_cast<std::size_t>(total_elems));
        DevicePtr<cuDoubleComplex> dK_ref(static_cast<std::size_t>(total_elems));
        if (opt.use_fused_f2k2) {
            v14_build_F2_batch_same_shape_kernel<<<blocks, threads, 0, stream>>>(
                dF2_ref.p, B, n, dim1, dim2, dEn.p,
                dp1x.p, dp1y.p, dp1z.p, de1l.p, de1m.p,
                dp2x.p, dp2y.p, dp2z.p, de2l.p, de2m.p,
                total_P[0].real(), total_P[1].real(), total_P[2].real(),
                par.atmK, par.atmpi, L, par.alpha, par.epsilon_h,
                (int)par.max_shell_num, par.Q0norm, 100, 1.0e-15);
            V11_CUDA_CHECK(cudaGetLastError());
            v14_build_K2_batch_same_shape_kernel<<<blocks, threads, 0, stream>>>(
                dK_ref.p, B, n, dim1, dim2, dEn.p,
                dp1x.p, dp1y.p, dp1z.p, de1l.p, de1m.p,
                dp2x.p, dp2y.p, dp2z.p, de2l.p, de2m.p,
                total_P[0].real(), total_P[1].real(), total_P[2].real(),
                par.eta_1, par.eta_2, dSP1.view(), dSP2.view(),
                par.atmK, par.atmpi, par.epsilon_h, L);
            V11_CUDA_CHECK(cudaGetLastError());
        } else {
            v32zu_fused_F2_K2inv_batch_same_shape_kernel<<<fused_blocks, threads, 0, stream>>>(
                dF2_ref.p, dK_ref.p, B, n, dim1, dim2, dEn.p,
                dp1x.p, dp1y.p, dp1z.p, de1l.p, de1m.p,
                dp2x.p, dp2y.p, dp2z.p, de2l.p, de2m.p,
                total_P[0].real(), total_P[1].real(), total_P[2].real(),
                par.eta_1, par.eta_2, dSP1.view(), dSP2.view(),
                par.atmK, par.atmpi, L, par.alpha, par.epsilon_h,
                (int)par.max_shell_num, par.Q0norm, 100, 1.0e-15);
            V11_CUDA_CHECK(cudaGetLastError());
        }
        sync_point();
        std::vector<cuDoubleComplex> hF2_ref(static_cast<std::size_t>(total_elems)), hF2_got(static_cast<std::size_t>(total_elems));
        std::vector<cuDoubleComplex> hK_ref(static_cast<std::size_t>(total_elems)), hK_got(static_cast<std::size_t>(total_elems));
        V11_CUDA_CHECK(cudaMemcpy(hF2_ref.data(), dF2_ref.p, sizeof(cuDoubleComplex)*hF2_ref.size(), cudaMemcpyDeviceToHost));
        V11_CUDA_CHECK(cudaMemcpy(hF2_got.data(), dF2.p, sizeof(cuDoubleComplex)*hF2_got.size(), cudaMemcpyDeviceToHost));
        V11_CUDA_CHECK(cudaMemcpy(hK_ref.data(), dK_ref.p, sizeof(cuDoubleComplex)*hK_ref.size(), cudaMemcpyDeviceToHost));
        V11_CUDA_CHECK(cudaMemcpy(hK_got.data(), dK.p, sizeof(cuDoubleComplex)*hK_got.size(), cudaMemcpyDeviceToHost));
        double aF=0.0,rF=0.0,aK=0.0,rK=0.0;
        v32zu_update_diff_stats(hF2_ref, hF2_got, aF, rF);
        v32zu_update_diff_stats(hK_ref, hK_got, aK, rK);
        g_v32zu_f2k2_validation.max_abs_diff_F2 = std::max(g_v32zu_f2k2_validation.max_abs_diff_F2, aF);
        g_v32zu_f2k2_validation.max_rel_diff_F2 = std::max(g_v32zu_f2k2_validation.max_rel_diff_F2, rF);
        g_v32zu_f2k2_validation.max_abs_diff_K2inv = std::max(g_v32zu_f2k2_validation.max_abs_diff_K2inv, aK);
        g_v32zu_f2k2_validation.max_rel_diff_K2inv = std::max(g_v32zu_f2k2_validation.max_rel_diff_K2inv, rK);
        g_v32zu_f2k2_validation.chunks += 1;
        std::cout << "[F2K2 validation] max_abs_diff_F2=" << aF
                  << " max_rel_diff_F2=" << rF
                  << " max_abs_diff_K2inv=" << aK
                  << " max_rel_diff_K2inv=" << rK << "\n";
    }

    t.reset();
    ggpu::Vec3d Pgg{total_P[0].real(), total_P[1].real(), total_P[2].real()};
    v14_build_G_batch_same_shape_kernel<<<blocks, threads, 0, stream>>>(
        dG.p, B, n, dim1, dim2, dEn.p, dg1.p, dg2.p, Pgg,
        par.atmK, par.atmpi, L, par.epsilon_h, par.Q0norm);
    V11_CUDA_CHECK(cudaGetLastError());
    v14_combine_H_batch_kernel<<<blocks, threads, 0, stream>>>(dH.p, dF2.p, dG.p, dK.p, total_elems);
    V11_CUDA_CHECK(cudaGetLastError());

    // H X1 = F2, build F3 = F2/3 - F2 X1, then solve F3 X = I.
    cublasHandle_t cublas = external_cublas;
    if (!cublas) V11_CUBLAS_CHECK(cublasCreate(&cublas));
    if (streamed_mode) CHECK_CUBLAS(cublasSetStream(cublas, stream));

    DevicePtr<cuDoubleComplex> dX1(static_cast<std::size_t>(total_elems));
    V11_CUDA_CHECK(cudaMemcpyAsync(dX1.p, dF2.p, sizeof(cuDoubleComplex)*static_cast<std::size_t>(total_elems), cudaMemcpyDeviceToDevice, stream));

    std::vector<cuDoubleComplex*> Hptr(static_cast<std::size_t>(B)), X1ptr(static_cast<std::size_t>(B));
    for (int e=0; e<B; ++e) {
        Hptr[static_cast<std::size_t>(e)]  = dH.p + static_cast<std::size_t>(e)*n*n;
        X1ptr[static_cast<std::size_t>(e)] = dX1.p   + static_cast<std::size_t>(e)*n*n;
    }

    sync_point();
    double gf3_part1_sec = t.elapsed();
    // We count G/H build in G/F3 matrix build time; the solve below is counted separately.

    t.reset();
    if (streamed_mode) v32zu_cublas_lu_solve_inplace_batched_same_n_stream(cublas, stream, n, n, Hptr, X1ptr, "H", debug);
    else cublas_lu_solve_inplace_batched_same_n(cublas, n, n, Hptr, X1ptr, debug);
    double solve_h_sec = t.elapsed();

    t.reset();
    DevicePtr<cuDoubleComplex*> dF2ptr(static_cast<std::size_t>(B)), dX1ptr(static_cast<std::size_t>(B)), dF2X1ptr(static_cast<std::size_t>(B));
    std::vector<cuDoubleComplex*> F2ptr(static_cast<std::size_t>(B)), F2X1ptr(static_cast<std::size_t>(B));
    for (int e=0; e<B; ++e) {
        F2ptr[static_cast<std::size_t>(e)]   = dF2.p   + static_cast<std::size_t>(e)*n*n;
        F2X1ptr[static_cast<std::size_t>(e)] = dH.p    + static_cast<std::size_t>(e)*n*n;
    }
    V11_CUDA_CHECK(cudaMemcpyAsync(dF2ptr.p, F2ptr.data(), sizeof(cuDoubleComplex*)*B, cudaMemcpyHostToDevice, stream));
    V11_CUDA_CHECK(cudaMemcpyAsync(dX1ptr.p, X1ptr.data(), sizeof(cuDoubleComplex*)*B, cudaMemcpyHostToDevice, stream));
    V11_CUDA_CHECK(cudaMemcpyAsync(dF2X1ptr.p, F2X1ptr.data(), sizeof(cuDoubleComplex*)*B, cudaMemcpyHostToDevice, stream));
    const cuDoubleComplex alpha1 = make_cuDoubleComplex(1.0,0.0);
    const cuDoubleComplex beta0  = make_cuDoubleComplex(0.0,0.0);
    V11_CUBLAS_CHECK(cublasZgemmBatched(cublas, CUBLAS_OP_N, CUBLAS_OP_N,
                                        n, n, n, &alpha1,
                                        const_cast<const cuDoubleComplex**>(dF2ptr.p), n,
                                        const_cast<const cuDoubleComplex**>(dX1ptr.p), n,
                                        &beta0, dF2X1ptr.p, n, B));
    v14_build_F3_batch_kernel<<<blocks, threads, 0, stream>>>(dH.p, dF2.p, dH.p, total_elems);
    V11_CUDA_CHECK(cudaGetLastError());
    sync_point();
    double gf3_part2_sec = t.elapsed();
    timing.gf3_build_sec += (gf3_part1_sec + gf3_part2_sec);
    std::cout << "[timing chunk] G/F3 matrix build time = " << (gf3_part1_sec + gf3_part2_sec) << " sec\n";

    t.reset();
    DevicePtr<cuDoubleComplex> dF3inv(static_cast<std::size_t>(total_elems));
    int iblocks = static_cast<int>((total_elems + threads - 1)/threads);
    v32zt_identity_batch_kernel<<<iblocks, threads, 0, stream>>>(dF3inv.p, B, n);
    V11_CUDA_CHECK(cudaGetLastError());

    std::vector<cuDoubleComplex*> F3ptr(static_cast<std::size_t>(B)), Invptr(static_cast<std::size_t>(B));
    for (int e=0; e<B; ++e) {
        F3ptr[static_cast<std::size_t>(e)] = dH.p      + static_cast<std::size_t>(e)*n*n;
        Invptr[static_cast<std::size_t>(e)] = dF3inv.p + static_cast<std::size_t>(e)*n*n;
    }
    if (streamed_mode) v32zu_cublas_lu_solve_inplace_batched_same_n_stream(cublas, stream, n, n, F3ptr, Invptr, "F3inv", debug);
    else cublas_lu_solve_inplace_batched_same_n(cublas, n, n, F3ptr, Invptr, debug);
    sync_point();
    double solve_f3_sec = t.elapsed();
    timing.solve_sec += (solve_h_sec + solve_f3_sec);
    std::cout << "[timing chunk] matrix inverse/solve time = " << (solve_h_sec + solve_f3_sec)
              << " sec (H solve=" << solve_h_sec << ", F3inv solve=" << solve_f3_sec << ")\n";

    // ---------------- Copy F3inv and Vsel to host, write records ----------------
    t.reset();
    std::vector<cuDoubleComplex> hF2(static_cast<std::size_t>(B)*n*n);
    std::vector<cuDoubleComplex> hK(static_cast<std::size_t>(B)*n*n);
    std::vector<cuDoubleComplex> hG(static_cast<std::size_t>(B)*n*n);
    std::vector<cuDoubleComplex> hF3(static_cast<std::size_t>(B)*n*n);
    std::vector<cuDoubleComplex> hF3inv(static_cast<std::size_t>(B)*n*n);
    std::vector<cuDoubleComplex> hVsel(static_cast<std::size_t>(B)*n*vdim);
    V11_CUDA_CHECK(cudaMemcpyAsync(hF2.data(), dF2.p, sizeof(cuDoubleComplex)*hF2.size(), cudaMemcpyDeviceToHost, stream));
    V11_CUDA_CHECK(cudaMemcpyAsync(hK.data(), dK.p, sizeof(cuDoubleComplex)*hK.size(), cudaMemcpyDeviceToHost, stream));
    V11_CUDA_CHECK(cudaMemcpyAsync(hG.data(), dG.p, sizeof(cuDoubleComplex)*hG.size(), cudaMemcpyDeviceToHost, stream));
    V11_CUDA_CHECK(cudaMemcpyAsync(hF3.data(), dH.p, sizeof(cuDoubleComplex)*hF3.size(), cudaMemcpyDeviceToHost, stream));
    V11_CUDA_CHECK(cudaMemcpyAsync(hF3inv.data(), dF3inv.p, sizeof(cuDoubleComplex)*hF3inv.size(), cudaMemcpyDeviceToHost, stream));
    V11_CUDA_CHECK(cudaMemcpyAsync(hVsel.data(), dVselBatch.p, sizeof(cuDoubleComplex)*hVsel.size(), cudaMemcpyDeviceToHost, stream));
    sync_point();

    for (int e=0; e<B; ++e)
    {
        const std::size_t off = static_cast<std::size_t>(e) * static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
        const cuDoubleComplex* oneF2 = hF2.data() + off;
        const cuDoubleComplex* oneK = hK.data() + off;
        const cuDoubleComplex* oneG = hG.data() + off;
        const cuDoubleComplex* oneF3 = hF3.data() + off;
        const cuDoubleComplex* oneF = hF3inv.data() + static_cast<std::size_t>(e) * static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
        const cuDoubleComplex* oneV = hVsel.data() + static_cast<std::size_t>(e) * static_cast<std::size_t>(n) * static_cast<std::size_t>(vdim);
        if (staged_rows)
        {
            V32ztCacheRowPayload row;
            row.grid_i = group[static_cast<std::size_t>(e)].i;
            row.Ecm = group[static_cast<std::size_t>(e)].Ecm;
            row.En = group[static_cast<std::size_t>(e)].En;
            row.dim1 = dim1;
            row.dim2 = dim2;
            row.n = n;
            row.vdim = vdim;
            row.F3inv_colmajor.assign(oneF, oneF + static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
            row.Vsel_colmajor.assign(oneV, oneV + static_cast<std::size_t>(n) * static_cast<std::size_t>(vdim));
            staged_rows->push_back(std::move(row));
        }
        if (allobj)
        {
            const auto row_i = group[static_cast<std::size_t>(e)].i;
            const auto row_Ecm = group[static_cast<std::size_t>(e)].Ecm;

            // Component-cache mode: write the actual GPU-produced matrices copied above.
            // F2/K2inv/G/F3/F3inv/Vsel are not recomputed on the CPU here.  The only
            // host-side matrix product is projected_F3inv = Vsel^dagger F3inv Vsel,
            // which is a cache convenience derived from the GPU F3inv and GPU Vsel.
            const Eigen::MatrixXcd plmM = v32zt_config_to_matrix(plm_all[static_cast<std::size_t>(e)]);
            const Eigen::MatrixXcd klmM = v32zt_config_to_matrix(klm_all[static_cast<std::size_t>(e)]);
            const Eigen::MatrixXcd F2M = v32zt_complex_buffer_to_matrix(oneF2, n, n);
            const Eigen::MatrixXcd KM = v32zt_complex_buffer_to_matrix(oneK, n, n);
            const Eigen::MatrixXcd GM = v32zt_complex_buffer_to_matrix(oneG, n, n);
            const Eigen::MatrixXcd F3M = v32zt_complex_buffer_to_matrix(oneF3, n, n);
            const Eigen::MatrixXcd F3invM = v32zt_complex_buffer_to_matrix(oneF, n, n);
            const Eigen::MatrixXcd VselM = v32zt_complex_buffer_to_matrix(oneV, n, vdim);
            const Eigen::MatrixXcd projectedF3invM = VselM.adjoint() * F3invM * VselM;

            v32zt_write_matrix_record(allobj->get("plm_config").os, row_i, row_Ecm, plmM);
            v32zt_write_matrix_record(allobj->get("klm_config").os, row_i, row_Ecm, klmM);
            v32zt_write_matrix_record(allobj->get("K2inv").os, row_i, row_Ecm, KM);
            v32zt_write_matrix_record(allobj->get("F2").os, row_i, row_Ecm, F2M);
            v32zt_write_matrix_record(allobj->get("G").os, row_i, row_Ecm, GM);
            v32zt_write_matrix_record(allobj->get("F3").os, row_i, row_Ecm, F3M);
            v32zt_write_matrix_record(allobj->get("F3inv").os, row_i, row_Ecm, F3invM);
            v32zt_write_matrix_record(allobj->get("Vsel").os, row_i, row_Ecm, VselM);
            v32zt_write_matrix_record(allobj->get("projected_F3inv").os, row_i, row_Ecm, projectedF3invM);
        }
        else if (staged_allobj_rows)
        {
            V32ztAllObjRowPayload row;
            row.grid_i = group[static_cast<std::size_t>(e)].i;
            row.Ecm = group[static_cast<std::size_t>(e)].Ecm;
            row.plm_config = v32zt_config_to_matrix(plm_all[static_cast<std::size_t>(e)]);
            row.klm_config = v32zt_config_to_matrix(klm_all[static_cast<std::size_t>(e)]);
            row.K2inv = v32zt_complex_buffer_to_matrix(oneK, n, n);
            row.F2 = v32zt_complex_buffer_to_matrix(oneF2, n, n);
            row.G = v32zt_complex_buffer_to_matrix(oneG, n, n);
            row.F3 = v32zt_complex_buffer_to_matrix(oneF3, n, n);
            row.F3inv = v32zt_complex_buffer_to_matrix(oneF, n, n);
            row.Vsel = v32zt_complex_buffer_to_matrix(oneV, n, vdim);
            row.projected_F3inv = row.Vsel.adjoint() * row.F3inv * row.Vsel;
            staged_allobj_rows->push_back(std::move(row));
        }
    }
    if (!staged_rows) fout.flush();
    double write_sec = t.elapsed();
    timing.write_sec += write_sec;
    std::cout << "[timing chunk] binary write time = " << write_sec << " sec\n";

    if (!external_cublas) V11_CUBLAS_CHECK(cublasDestroy(cublas));

    timing.rows += static_cast<std::size_t>(B);
    timing.chunks += 1;
    timing.max_total_dim = std::max(timing.max_total_dim, n);
    const double pct = 100.0 * double(global_done_before + static_cast<std::size_t>(B)) / std::max<std::size_t>(1, global_total);
    std::cout << std::fixed << std::setprecision(1)
              << "[progress] irrep=" << spec.label
              << " rows_done=" << (global_done_before + static_cast<std::size_t>(B)) << "/" << global_total
              << " (" << pct << "%) chunk_time=" << t_total.elapsed() << " sec\n"
              << std::setprecision(6) << std::defaultfloat;
}

static void v32zu_write_staged_rows(
    std::ofstream& fout,
    const std::vector<V32ztCacheRowPayload>& rows,
    V32ztTiming& timing)
{
    // v32zu streams: serialize completed concurrent rows with the original binary record layout.
    V32ztStopwatch twrite;
    for (const auto& row : rows)
    {
        v32zt_write_cache_record(fout,
                                  row.grid_i,
                                  row.Ecm,
                                  row.En,
                                  row.dim1,
                                  row.dim2,
                                  row.n,
                                  row.vdim,
                                  row.F3inv_colmajor,
                                  row.Vsel_colmajor,
                                  timing.bytes_written);
    }
    fout.flush();
    timing.write_sec += twrite.elapsed();
}

static void v32zu_accumulate_timing_without_write(V32ztTiming& dst, const V32ztTiming& src)
{
    // v32zu streams: compute timings are produced per worker; binary write timing is measured during ordered flush.
    dst.basis_sec += src.basis_sec;
    dst.f2k2_sec += src.f2k2_sec;
    dst.gf3_build_sec += src.gf3_build_sec;
    dst.solve_sec += src.solve_sec;
    dst.vsel_sec += src.vsel_sec;
    dst.rows += src.rows;
    dst.chunks += src.chunks;
    dst.max_total_dim = std::max(dst.max_total_dim, src.max_total_dim);
    dst.max_vdim = std::max(dst.max_vdim, src.max_vdim);
}

static void run_groups_concurrent(
    const std::map<EnergyShapeKey, std::vector<EnergyPrepassInfo>>& groups,
    const V32ztIrrepSpec& spec,
    const PhysicsParams& par,
    std::ofstream& fout,
    V32ztTiming& timing,
    std::size_t global_total,
    int total_chunks_planned,
    const V32zuRunOptions& opt,
    char debug,
    V32ztAllObjWriter* allobj_writer)
{
    // v32zu streams: schedule same-size chunks on separate streams while preserving original write order.
    std::vector<V32zuChunkTask> tasks;
    int group_id = 0;
    int chunk_id = 0;
    int original_order = 0;
    for (const auto& kv : groups)
    {
        ++group_id;
        const int n = kv.first.total_dim;
        const auto& group = kv.second;
        const int cap = v32zt_auto_batch_cap(n, opt.max_batch_energies, opt.gpu_memory_safety_fraction);
        std::size_t free_b=0,total_b=0; cudaMemGetInfo(&free_b,&total_b);
        const std::size_t estimated_per_row = v32zt_estimate_group_bytes_per_energy(n);
        std::cout << "[VRAM] irrep=" << spec.label << " n=" << n
                  << " free_GB=" << (double(free_b)/1.0e9)
                  << " total_GB=" << (double(total_b)/1.0e9)
                  << " safety_fraction=" << opt.gpu_memory_safety_fraction
                  << " concurrent_streams_safety_fraction=" << opt.concurrent_streams_safety_fraction
                  << " estimated_bytes_per_row=" << estimated_per_row
                  << " max_batch_for_group=" << cap << "\n";
        for (std::size_t begin=0; begin<group.size(); begin += static_cast<std::size_t>(cap))
        {
            const std::size_t end = std::min(group.size(), begin + static_cast<std::size_t>(cap));
            V32zuChunkTask task;
            task.original_order = original_order++;
            task.group_id = group_id;
            task.chunk_id = ++chunk_id;
            task.dim1 = kv.first.dim1;
            task.dim2 = kv.first.dim2;
            task.n = n;
            task.rows.assign(group.begin()+static_cast<long>(begin), group.begin()+static_cast<long>(end));
            task.estimated_vram_bytes = estimate_group_vram_bytes(task.rows);
            tasks.push_back(std::move(task));
        }
    }

    std::vector<int> unfinished(tasks.size());
    for (std::size_t i=0; i<tasks.size(); ++i) unfinished[i] = static_cast<int>(i);
    std::map<int, V32zuChunkResult> completed_by_order;
    std::size_t rows_flushed = 0;
    int next_write_order = 0;
    int wave_id = 0;

    auto flush_ready = [&]() {
        while (true)
        {
            auto it = completed_by_order.find(next_write_order);
            if (it == completed_by_order.end()) break;
            v32zu_write_staged_rows(fout, it->second.rows, timing);
            if (allobj_writer && !it->second.allobj_rows.empty())
            {
                v32zt_write_staged_allobj_rows(*allobj_writer, it->second.allobj_rows, timing);
            }
            rows_flushed += it->second.rows.size();
            const double pct = 100.0 * double(rows_flushed) / std::max<std::size_t>(1, global_total);
            std::cout << std::fixed << std::setprecision(1)
                      << "[progress] irrep=" << spec.label
                      << " rows_done=" << rows_flushed << "/" << global_total
                      << " (" << pct << "%) ordered_flush_chunk=" << it->second.chunk_id << "\n"
                      << std::setprecision(6) << std::defaultfloat;
            completed_by_order.erase(it);
            ++next_write_order;
        }
    };

    while (!unfinished.empty())
    {
        std::sort(unfinished.begin(), unfinished.end(), [&](int a, int b) {
            return tasks[static_cast<std::size_t>(a)].estimated_vram_bytes > tasks[static_cast<std::size_t>(b)].estimated_vram_bytes;
        });

        const std::size_t free_now = get_free_vram_bytes();
        const std::size_t rss_now = v32zt_process_rss_bytes();
        const std::size_t swap_now = v32zt_process_swap_bytes();
        const std::size_t host_pressure_now = rss_now + swap_now;
        const std::size_t budget = static_cast<std::size_t>(double(free_now) * std::min(0.95, std::max(0.10, opt.concurrent_streams_safety_fraction)));
        std::vector<int> gpu_candidates;
        std::vector<int> serial_candidates;
        for (int idx : unfinished)
        {
            if (v32zt_should_use_gpu_lane(tasks[static_cast<std::size_t>(idx)], opt, free_now, host_pressure_now))
                gpu_candidates.push_back(idx);
            else
                serial_candidates.push_back(idx);
        }

        std::cout << "[scheduler] wave=" << (wave_id + 1)
                  << " free_GB=" << (double(free_now) / 1.0e9)
                  << " rss_GB=" << (double(rss_now) / 1.0e9)
                  << " swap_GB=" << (double(swap_now) / 1.0e9)
                  << " host_pressure_GB=" << (double(host_pressure_now) / 1.0e9)
                  << " gpu_candidates=" << gpu_candidates.size()
                  << " serial_candidates=" << serial_candidates.size()
                  << " budget_GB=" << (double(budget) / 1.0e9) << "\n";

        if (gpu_candidates.empty())
        {
            const int idx = unfinished.back();
            unfinished.pop_back();
            const V32zuChunkTask& task = tasks[static_cast<std::size_t>(idx)];
            std::cout << "[scheduler lane] chunk=" << task.chunk_id
                      << " lane=serial reason=insufficient_headroom_or_small_task\n";
            V32zuChunkResult result;
            result.original_order = task.original_order;
            result.chunk_id = task.chunk_id;
            result.n = task.n;
            result.dim1 = task.dim1;
            result.dim2 = task.dim2;
            V32ztTiming local_timing;
            std::ofstream dummy;
            v32zt_process_same_shape_chunk_cache(task.rows, spec, par, dummy, local_timing,
                                                 task.chunk_id, total_chunks_planned, 0, global_total,
                                                 opt, debug, &result.rows,
                                                 allobj_writer ? &result.allobj_rows : nullptr,
                                                 nullptr, nullptr);
            result.timing = local_timing;
            v32zu_accumulate_timing_without_write(timing, result.timing);
            completed_by_order[result.original_order] = std::move(result);
            flush_ready();
            continue;
        }

        std::vector<int> wave;
        std::size_t wave_bytes = 0;
        std::vector<int> rest;
        const std::size_t host_budget = v32zt_host_rss_budget_bytes();
        const std::size_t host_soft_limit = host_budget > 0 ? static_cast<std::size_t>(double(host_budget) * 0.80) : 0;
        const bool force_single_task_wave = allobj_writer != nullptr || (host_soft_limit > 0 && host_pressure_now >= host_soft_limit);
        if (force_single_task_wave)
        {
            // ponytail: one chunk at a time is the ceiling once writer staging or host pressure gets tight.
            const int idx = gpu_candidates.front();
            wave.push_back(idx);
            wave_bytes += tasks[static_cast<std::size_t>(idx)].estimated_vram_bytes;
            rest.insert(rest.end(), gpu_candidates.begin() + 1, gpu_candidates.end());
            std::cout << "[scheduler] wave=" << (wave_id + 1)
                      << " mode=single_task"
                      << " reason=" << (allobj_writer != nullptr ? "writer_bound" : "host_pressure")
                      << " host_soft_limit_GB=" << (double(host_soft_limit) / 1.0e9) << "\n";
        }
        else
        {
            for (int idx : gpu_candidates)
            {
                const std::size_t est = tasks[static_cast<std::size_t>(idx)].estimated_vram_bytes;
                if (wave.empty() || wave_bytes + est <= budget)
                {
                    wave.push_back(idx);
                    wave_bytes += est;
                }
                else
                {
                    rest.push_back(idx);
                }
            }
        }
        rest.insert(rest.end(), serial_candidates.begin(), serial_candidates.end());
        unfinished.swap(rest);
        ++wave_id;
        std::cout << "[wave " << wave_id << "] launching " << wave.size()
                  << " groups, combined VRAM = " << (double(wave_bytes)/1.0e9)
                  << " GB / " << (double(free_now)/1.0e9)
                  << " GB free, safety_fraction=" << opt.concurrent_streams_safety_fraction << "\n";

        std::vector<GroupStream> streams(wave.size());
        std::vector<std::shared_ptr<std::promise<void>>> ready_promises(wave.size());
        std::vector<std::future<void>> ready_futures;
        std::promise<void> launch_promise;
        std::shared_future<void> launch_gate(launch_promise.get_future());
        std::vector<std::future<V32zuChunkResult>> futures;

        for (std::size_t wi=0; wi<wave.size(); ++wi)
        {
            const V32zuChunkTask& task = tasks[static_cast<std::size_t>(wave[wi])];
            GroupStream& gs = streams[wi];
            gs.group_id = task.group_id;
            gs.n = task.n;
            gs.dim1 = task.dim1;
            gs.dim2 = task.dim2;
            gs.batch_size = static_cast<int>(task.rows.size());
            gs.estimated_vram_bytes = task.estimated_vram_bytes;
            CHECK_CUDA(cudaStreamCreateWithFlags(&gs.stream, cudaStreamNonBlocking));
            CHECK_CUBLAS(cublasCreate(&gs.handle));
            CHECK_CUBLAS(cublasSetStream(gs.handle, gs.stream));
            ready_promises[wi] = std::make_shared<std::promise<void>>();
            ready_futures.push_back(ready_promises[wi]->get_future());
            std::cout << "[wave group] wave=" << wave_id
                      << " group_id=" << gs.group_id
                      << " chunk=" << task.chunk_id
                      << " dim1=" << gs.dim1 << " dim2=" << gs.dim2
                      << " n=" << gs.n
                      << " batch_rows=" << gs.batch_size
                      << " estimated_GB=" << (double(gs.estimated_vram_bytes)/1.0e9)
                      << "\n";

            futures.push_back(std::async(std::launch::async,
                [&, task, wi]() mutable -> V32zuChunkResult {
                    try
                    {
                        V32zuChunkResult result;
                        result.original_order = task.original_order;
                        result.chunk_id = task.chunk_id;
                        result.n = task.n;
                        result.dim1 = task.dim1;
                        result.dim2 = task.dim2;
                        V32ztTiming local_timing;
                        std::ofstream dummy;
                        v32zt_process_same_shape_chunk_cache(task.rows, spec, par, dummy, local_timing,
                                                             task.chunk_id, total_chunks_planned, 0, global_total,
                                                             opt, debug, &result.rows,
                                                             allobj_writer ? &result.allobj_rows : nullptr,
                                                             nullptr,
                                                             streams[wi].stream, streams[wi].handle,
                                                             launch_gate, ready_promises[wi].get());
                        result.timing = local_timing;
                        streams[wi].completed = true;
                        return result;
                    }
                    catch (...)
                    {
                        try { ready_promises[wi]->set_exception(std::current_exception()); } catch (...) {}
                        throw;
                    }
                }));
            streams[wi].launched = true;
        }

        bool ready_ok = true;
        for (auto& rf : ready_futures)
        {
            try { rf.get(); }
            catch (const std::exception& ex)
            {
                ready_ok = false;
                std::cout << "[VRAM fallback] wave=" << wave_id
                          << " reason=prepare_failed detail=" << ex.what() << "\n";
            }
        }
        launch_promise.set_value();

        std::vector<int> failed_task_indices;
        for (std::size_t wi=0; wi<futures.size(); ++wi)
        {
            try
            {
                V32zuChunkResult result = futures[wi].get();
                v32zu_accumulate_timing_without_write(timing, result.timing);
                completed_by_order[result.original_order] = std::move(result);
            }
            catch (const std::exception& ex)
            {
                const V32zuChunkTask& task = tasks[static_cast<std::size_t>(wave[wi])];
                std::size_t free_b=0,total_b=0; cudaMemGetInfo(&free_b,&total_b);
                std::cout << "[VRAM fallback] wave=" << wave_id
                          << " n=" << task.n
                          << " dim1=" << task.dim1
                          << " dim2=" << task.dim2
                          << " batch_rows=" << task.rows.size()
                          << " estimated_bytes=" << task.estimated_vram_bytes
                          << " free_GB=" << (double(free_b)/1.0e9)
                          << " rss_GB=" << (double(v32zt_process_rss_bytes())/1.0e9)
                          << " safety_fraction=" << opt.concurrent_streams_safety_fraction
                          << " detail=" << ex.what()
                          << " action=retry_sequential\n";
                failed_task_indices.push_back(wave[wi]);
            }
        }
        flush_ready();

        for (GroupStream& gs : streams)
        {
            if (gs.stream) cudaStreamSynchronize(gs.stream);
            if (gs.handle) cublasDestroy(gs.handle);
            if (gs.stream) cudaStreamDestroy(gs.stream);
        }

        for (int failed_idx : failed_task_indices)
        {
            const V32zuChunkTask& task = tasks[static_cast<std::size_t>(failed_idx)];
            V32zuChunkResult result;
            result.original_order = task.original_order;
            result.chunk_id = task.chunk_id;
            result.n = task.n;
            result.dim1 = task.dim1;
            result.dim2 = task.dim2;
            V32ztTiming local_timing;
            std::ofstream dummy;
            v32zt_process_same_shape_chunk_cache(task.rows, spec, par, dummy, local_timing,
                                                 task.chunk_id, total_chunks_planned, 0, global_total,
                                                 opt, debug, &result.rows,
                                                 allobj_writer ? &result.allobj_rows : nullptr,
                                                 nullptr);
            result.timing = local_timing;
            v32zu_accumulate_timing_without_write(timing, result.timing);
            completed_by_order[result.original_order] = std::move(result);
        }
        flush_ready();
        if (!ready_ok) std::cout << "[wave " << wave_id << "] completed with prepare fallback diagnostics\n";
    }
    flush_ready();
    if (!completed_by_order.empty())
        throw std::runtime_error("concurrent scheduler internal error: unflushed staged chunks remain");
}

static V32ztTiming v32zt_generate_one_irrep(
    const V32ztIrrepSpec& spec,
    PhysicsParams par,
    double E0,
    double E1,
    int coarseN,
    const std::string& outdir,
    const V32zuRunOptions& opt,
    int assigned_device,
    char debug)
{
    if (assigned_device >= 0) {                                                          // v32zu OPT-8: sequential per-irrep device assignment keeps filenames/write ownership unchanged.
        V11_CUDA_CHECK(cudaSetDevice(assigned_device));
    }
#ifdef V32ZU_ALLOBJ_QC_WRITER
    const bool allobj_mode = true;
#else
    const bool allobj_mode = false;
#endif
    V32ztTiming timing;
    V32ztStopwatch irrep_timer;

    const std::string cache_name = std::string("v32zu_Lbyas") + v32zt_ltag(par.Lbyas) +
        "_xi" + v32zt_ltag(par.xi) + "_irrep" + spec.label + "_coarse" + std::to_string(coarseN) +
        "_F3inv_Vsel_gpu.bin";
    const std::string cache_path = v32zt_join_path(outdir, cache_name);
    const std::string meta_path = cache_path + ".meta.json";

    std::ofstream fout(cache_path, std::ios::binary);
    if (!fout) throw std::runtime_error("could not open cache output file: " + cache_path);
    std::size_t header_bytes = 0;
    const std::string header = "KKPI_F3INV_VSEL_GPU_V32ZT_EXACT_V32Y 1\n";
    v32zt_write_exact(fout, header.data(), header.size(), header_bytes);
    timing.bytes_written += header_bytes;

    std::cout << "============================================================\n";
    std::cout << "[IRREP START] " << spec.label << " Lbyas=" << par.Lbyas
              << " coarseN=" << coarseN << " Ecm=[" << std::setprecision(17) << E0 << "," << E1 << "]\n";
    std::cout << "[CACHE] " << cache_path << "\n";
    std::cout << "[pipeline] requested=" << (opt.use_stream_pipeline ? "true" : "false")
              << " enabled=" << (opt.use_concurrent_streams ? "true" : "false")
              << " reason=" << (opt.use_concurrent_streams ? "concurrent_same_size_group_waves" : "sequential_fallback_requested") << "\n";

    V32ztStopwatch tpre;
    const double del = (coarseN > 1) ? (E1 - E0) / double(coarseN - 1) : 0.0;
    std::vector<EnergyPrepassInfo> prepass_all(static_cast<std::size_t>(coarseN));
    int last_percent = -10;
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<coarseN; ++i)
    {
        const double Ecm = E0 + double(i) * del;
        prepass_all[static_cast<std::size_t>(i)] = v32zt_energy_shape_prepass_momentum_first(i, Ecm, spec.nnP, par, debug);
        int pct = (100 * (i + 1)) / std::max(1, coarseN);
        #pragma omp critical(v32zt_prepass_print)
        {
            if (pct >= last_percent + 10) {
                last_percent = pct;
                std::cout << "[prepass progress] irrep=" << spec.label << " " << pct << "%\n";
            }
        }
    }
    double prepass_sec = tpre.elapsed();
    std::cout << "[timing irrep] prepass/basis-shape construction time = " << prepass_sec << " sec\n";
    timing.basis_sec += prepass_sec;

    V32ztStopwatch tkin;                                                                  // v32zu OPT-4: only row energy/shape precompute is currently safe and already produced by prepass.
    const std::size_t kin_entries = prepass_all.size();
    const std::size_t kin_bytes = kin_entries * (sizeof(double) * 2 + sizeof(int) * 4);
    std::cout << "[kinematic table] enabled=" << (opt.use_kinematic_precompute ? "partial_safe" : "false")
              << " bytes=" << kin_bytes
              << " entries=" << kin_entries
              << " note=Ecm/En/shape prepass only; zeta/F2 helpers remain direct-evaluated\n";
    std::cout << "[timing irrep] kinematic precompute time = " << tkin.elapsed() << " sec\n";

    std::vector<EnergyPrepassInfo> prepass;
    prepass.reserve(prepass_all.size());
    for (const auto& info : prepass_all) {
        if (info.total_dim > 0) prepass.push_back(info);
        else std::cout << "[skip] empty config row i=" << info.i << " Ecm=" << info.Ecm << "\n";
    }

    std::unique_ptr<V32ztAllObjWriter> allobj_writer;
    if (allobj_mode) {
        allobj_writer = std::make_unique<V32ztAllObjWriter>(
            v32zt_open_allobj_writers(outdir, spec, par, coarseN, prepass.size()));
    }

    std::map<EnergyShapeKey, std::vector<EnergyPrepassInfo>> groups;
    for (const auto& info : prepass)
    {
        EnergyShapeKey key;
        key.dim1 = info.dim1;
        key.dim2 = info.dim2;
        key.total_dim = info.total_dim;
        key.vdim = -1;
        groups[key].push_back(info);
    }

    std::cout << "[groups] irrep=" << spec.label << " groups=" << groups.size()
              << " nonempty_rows=" << prepass.size() << "\n";
    for (const auto& kv : groups)
    {
        const auto& key = kv.first;
        const auto& group = kv.second;
        std::cout << "[group] dim1=" << key.dim1 << " dim2=" << key.dim2
                  << " n=" << key.total_dim << " count=" << group.size() << "\n";
    }

    std::size_t done = 0;
    int total_chunks_planned = 0;
    for (const auto& kv : groups) {
        int cap = v32zt_auto_batch_cap(kv.first.total_dim, opt.max_batch_energies, opt.gpu_memory_safety_fraction);
        total_chunks_planned += static_cast<int>((kv.second.size() + std::size_t(cap) - 1) / std::size_t(cap));
    }
    int chunk_id = 0;
#if USE_CONCURRENT_STREAMS
    if (opt.use_concurrent_streams)
    {
        run_groups_concurrent(groups, spec, par, fout, timing, prepass.size(), total_chunks_planned, opt, debug,
                              allobj_writer.get());
    }
    else
#endif
    for (const auto& kv : groups)
    {
        const int n = kv.first.total_dim;
        const auto& group = kv.second;
        const int cap = v32zt_auto_batch_cap(n, opt.max_batch_energies, opt.gpu_memory_safety_fraction);
        std::size_t free_b=0,total_b=0; cudaMemGetInfo(&free_b,&total_b);
        const std::size_t estimated_per_row = v32zt_estimate_group_bytes_per_energy(n);
        std::cout << "[VRAM] irrep=" << spec.label << " n=" << n
                  << " free_GB=" << (double(free_b)/1.0e9)
                  << " total_GB=" << (double(total_b)/1.0e9)
                  << " safety_fraction=" << opt.gpu_memory_safety_fraction
                  << " estimated_bytes_per_row=" << estimated_per_row
                  << " max_batch_for_group=" << cap << "\n";
        for (std::size_t begin=0; begin<group.size(); )
        {
            std::size_t current_cap = static_cast<std::size_t>(cap);                      // v32zu OPT-6: adapt this group chunk if cudaMalloc fails.
            const std::size_t end = std::min(group.size(), begin + current_cap);
            std::vector<EnergyPrepassInfo> chunk(group.begin()+static_cast<long>(begin), group.begin()+static_cast<long>(end));
            while (!chunk.empty()) {
                chunk_id++;
                try {
                    v32zt_process_same_shape_chunk_cache(chunk, spec, par, fout, timing,
                                                         chunk_id, total_chunks_planned, done, prepass.size(), opt, debug,
                                                         nullptr, nullptr, allobj_writer.get());
                    done += chunk.size();
                    begin += chunk.size();
                    break;
                } catch (const std::exception& ex) {
                    const std::string msg = ex.what();
                    const bool malloc_failed = (msg.find("cudaMalloc") != std::string::npos) || (msg.find("out of memory") != std::string::npos);
                    std::cout << "[VRAM fallback] n=" << n
                              << " dim1=" << kv.first.dim1
                              << " dim2=" << kv.first.dim2
                              << " requested_batch_rows=" << chunk.size()
                              << " estimated_bytes=" << (estimated_per_row * chunk.size())
                              << " free_GB=" << (double(free_b)/1.0e9)
                              << " reason=" << (malloc_failed ? "cudaMalloc_failed" : "chunk_failed")
                              << " detail=" << msg << "\n";
                    cudaDeviceSynchronize();                                              // v32zu OPT-6: clear device state before retrying smaller chunks without aborting fallback.
                    cudaGetLastError();
                    if (chunk.size() == 1) throw;
                    const std::size_t old_batch = chunk.size();
                    const std::size_t new_batch = std::max<std::size_t>(1, old_batch / 2);
                    std::cout << "[VRAM fallback] old_batch=" << old_batch
                              << " new_batch=" << new_batch
                              << " reason=" << (malloc_failed ? "cudaMalloc_failed" : "chunk_failed") << "\n";
                    chunk.assign(group.begin()+static_cast<long>(begin),
                                 group.begin()+static_cast<long>(begin + new_batch));
                }
            }
        }
    }

    timing.total_sec = irrep_timer.elapsed();

    if (allobj_writer) {
        v32zt_close_allobj_writers(*allobj_writer);
        cudaDeviceProp prop{};
        int dev = 0;
        if (cudaGetDevice(&dev) != cudaSuccess) dev = 0;
        if (cudaGetDeviceProperties(&prop, dev) != cudaSuccess) {
            std::memset(&prop, 0, sizeof(prop));
            std::snprintf(prop.name, sizeof(prop.name), "%s", "unknown");
        }
        v32zt_finalize_allobj_sidecars(*allobj_writer, spec, par, E0, E1, coarseN, timing, dev, prop,
                                       opt.use_concurrent_streams);
    }

    fout.close();
    v32zt_write_sidecar_json(meta_path, spec, par, E0, E1, coarseN, cache_path, timing, opt, assigned_device,
                             allobj_writer != nullptr && opt.use_concurrent_streams);

    std::cout << "============================================================\n";
    std::cout << "[IRREP SUMMARY] " << spec.label << "\n";
    std::cout << "rows written                  = " << timing.rows << "\n";
    std::cout << "chunks                        = " << timing.chunks << "\n";
    std::cout << "max total dim                 = " << timing.max_total_dim << "\n";
    std::cout << "max vdim                      = " << timing.max_vdim << "\n";
    std::cout << "basis construction time       = " << timing.basis_sec << " sec\n";
    std::cout << "F2/K2 time                    = " << timing.f2k2_sec << " sec\n";
    std::cout << "G/F3 matrix build time        = " << timing.gf3_build_sec << " sec\n";
    std::cout << "matrix inverse/solve time     = " << timing.solve_sec << " sec\n";
    std::cout << "Vsel construction time        = " << timing.vsel_sec << " sec\n";
    std::cout << "binary write time             = " << timing.write_sec << " sec\n";
    std::cout << "total irrep time              = " << timing.total_sec << " sec (" << v32zt_format_seconds(timing.total_sec) << ")\n";
    std::cout << "cache                         = " << cache_path << "\n";
    std::cout << "sidecar                       = " << meta_path << "\n";
    return timing;
}

static void v32zt_print_usage(const char* prog)
{
    std::cerr << "Usage:\n"
              << "  " << prog << " <Lbyas> <xi> <Ecm_min> <Ecm_max> <coarseN> <output_dir> [debug] [max_batch_energies] [gpu_memory_safety_fraction] [use_fused_f2k2] [validate_fused_f2k2] [use_kinematic_precompute] [use_stream_pipeline] [use_multi_gpu] [use_mixed_precision_f2k2] [use_concurrent_streams] [concurrent_streams_safety_fraction] [irreps_csv]\n";
}

int main(int argc, char** argv)
{
    std::signal(SIGSEGV, v11_sigsegv_handler);
    try {
        if (argc < 7) { v32zt_print_usage(argv[0]); return 2; }

        PhysicsParams par;
        // Strict v32y stream-cache convention: global parity=-1 for all five default irreps.
        par.parity = -1;
        par.Lbyas = std::stod(argv[1]);
        par.xi = std::stod(argv[2]);
        const double E0 = std::stod(argv[3]);
        const double E1 = std::stod(argv[4]);
        const int coarseN = std::stoi(argv[5]);
        const std::string outdir = argv[6];
        const char debug = (argc >= 8 && !std::string(argv[7]).empty()) ? argv[7][0] : 'n';
        V32zuRunOptions opt;                                                             // v32zu OPT-1..8: parse optional flags after the legacy positional args.
        opt.max_batch_energies = (argc >= 9) ? std::stoi(argv[8]) : 0;
        opt.gpu_memory_safety_fraction = (argc >= 10) ? std::stod(argv[9]) : 0.85;
        opt.use_fused_f2k2 = (argc >= 11) ? v32zu_parse_bool(argv[10]) : true;
        opt.validate_fused_f2k2 = (argc >= 12) ? v32zu_parse_bool(argv[11]) : false;
        opt.use_kinematic_precompute = (argc >= 13) ? v32zu_parse_bool(argv[12]) : true;
        opt.use_stream_pipeline = (argc >= 14) ? v32zu_parse_bool(argv[13]) : true;
        opt.use_multi_gpu = (argc >= 15) ? v32zu_parse_bool(argv[14]) : false;
        opt.use_mixed_precision_f2k2 = (argc >= 16) ? v32zu_parse_bool(argv[15]) : false;
        opt.use_concurrent_streams = (argc >= 17) ? v32zu_parse_bool(argv[16]) : true;
        opt.concurrent_streams_safety_fraction = (argc >= 18) ? std::stod(argv[17]) : CONCURRENT_STREAMS_VRAM_SAFETY_FRACTION;
        opt.irrep_filter_csv = (argc >= 19) ? std::string(argv[18]) : "";
#ifdef V32ZU_ALLOBJ_QC_WRITER
        opt.use_concurrent_streams = true;
#endif
#if !USE_CONCURRENT_STREAMS
        if (opt.use_concurrent_streams)
            std::cout << "[streams] requested concurrent streams but binary was compiled with USE_CONCURRENT_STREAMS=0; using sequential fallback\n";
        opt.use_concurrent_streams = false;
#endif
        opt.concurrent_streams_safety_fraction = std::min(0.95, std::max(0.10, opt.concurrent_streams_safety_fraction));
#ifdef USE_MIXED_PRECISION_F2K2
        const bool mixed_precision_compiled = true;                                      // v32zu OPT-7: compile-time guard remains off unless explicitly defined.
#else
        const bool mixed_precision_compiled = false;
#endif
        if (opt.use_mixed_precision_f2k2 && !mixed_precision_compiled) {
            std::cout << "[precision] requested mixed precision but binary was not compiled with USE_MIXED_PRECISION_F2K2; using FP64\n";
            opt.use_mixed_precision_f2k2 = false;
        }

        // Hard safety check requested by Digonto.
        const double scatter1_10 = -43.2;
        if (!std::isfinite(scatter1_10) || std::abs(scatter1_10) < 1.0e-300)
            throw std::runtime_error("p-wave ell=1 is active, but scatter1_10 is zero/non-finite");

        V11_CUDA_CHECK(cudaGetDeviceCount(&opt.device_count));                           // v32zu OPT-8: detect devices up front for guarded multi-GPU assignment.
        int dev = 0;
        V11_CUDA_CHECK(cudaGetDevice(&dev));
        cudaDeviceProp prop{};
        V11_CUDA_CHECK(cudaGetDeviceProperties(&prop, dev));
        std::size_t free_b=0,total_b=0; cudaMemGetInfo(&free_b,&total_b);

        std::cout << "# ============================================================\n";
        std::cout << "# " << V32ZT_VERSION << "\n";
        std::cout << "# Lbyas=" << std::setprecision(17) << par.Lbyas
                  << " xi=" << par.xi << " spatial_L=" << par.xi*par.Lbyas << "\n";
        std::cout << "# Ecm_min=" << E0 << " Ecm_max=" << E1 << " coarseN=" << coarseN << "\n";
        std::cout << "# output_dir=" << outdir << "\n";
        std::cout << "# default irreps: 000_A1m 100_A2 110_A2 111_A2 200_A2\n";
        std::cout << "# scatter safe defaults: scatter1_00=4.04 scatter1_10=-43.2 scatter2_00=4.12\n";
        std::cout << "# waves_vec_1={0,1}, waves_vec_2={0}\n";
        std::cout << "# Vsel projector=canonical k3iso01 CPU convention; parity=-1 eig_tol=0.05\n";
        std::cout << "[config] use_fused_f2k2=" << (opt.use_fused_f2k2 ? "true" : "false")
                  << " validate_fused_f2k2=" << (opt.validate_fused_f2k2 ? "true" : "false")
                  << " use_kinematic_precompute=" << (opt.use_kinematic_precompute ? "true" : "false")
                  << " use_stream_pipeline=" << (opt.use_stream_pipeline ? "true" : "false")
                  << " use_concurrent_streams=" << (opt.use_concurrent_streams ? "true" : "false")
                  << " max_batch_energies=" << opt.max_batch_energies
                  << " gpu_memory_safety_fraction=" << opt.gpu_memory_safety_fraction
                  << " concurrent_streams_safety_fraction=" << opt.concurrent_streams_safety_fraction
                  << " use_multi_gpu=" << (opt.use_multi_gpu ? "true" : "false") << "\n";
        std::cout << "[precision] F2/K2 mixed precision enabled=" << (opt.use_mixed_precision_f2k2 ? "true" : "false")
                  << " compiled_flag=" << (mixed_precision_compiled ? "true" : "false")
                  << " default_path=FP64\n";
        std::cout << "[multiGPU] device_count=" << opt.device_count
                  << " requested=" << (opt.use_multi_gpu ? "true" : "false");
        if (!opt.use_multi_gpu) std::cout << " enabled=false reason=config_disabled\n";
        else if (opt.device_count <= 1) std::cout << " enabled=false only one GPU found; using single GPU\n";
        else std::cout << " enabled=true mode=sequential_irrep_assignment\n";
        std::cout << "[streams] concurrent_streams enabled=" << (opt.use_concurrent_streams ? "true" : "false") << "\n";
        std::cout << "[streams] safety_fraction=" << opt.concurrent_streams_safety_fraction << "\n";
        std::cout << "[streams] fallback sequential path available=true\n";
        std::cout << "# GPU device=" << dev << " name=" << prop.name
                  << " free_GB=" << (double(free_b)/1.0e9)
                  << " total_GB=" << (double(total_b)/1.0e9) << "\n";
        std::cout << "# max_batch_energies=" << opt.max_batch_energies << " (0 means auto), safety_fraction=" << opt.gpu_memory_safety_fraction << "\n";
        std::cout << "[F2K2 staging] mode=register_read_only table_bytes=0 note=row energy and constants staged in registers; larger tables use global fallback\n";
        std::cout << "# ============================================================\n";

        V32ztStopwatch Ltimer;
        std::vector<V32ztTiming> irrep_timings;
        auto irreps = v32zt_default_irreps();
        if (!opt.irrep_filter_csv.empty())
        {
            // v32zu streams: smoke tests can restrict irreps without changing the default production set.
            std::vector<std::string> wanted;
            std::stringstream ss(opt.irrep_filter_csv);
            std::string item;
            while (std::getline(ss, item, ',')) if (!item.empty()) wanted.push_back(item);
            std::vector<V32ztIrrepSpec> filtered;
            for (const auto& ir : irreps)
                if (std::find(wanted.begin(), wanted.end(), ir.label) != wanted.end()) filtered.push_back(ir);
            if (filtered.empty())
                throw std::runtime_error("irrep filter selected no known irreps: " + opt.irrep_filter_csv);
            irreps = filtered;
            std::cout << "[config] irrep_filter=" << opt.irrep_filter_csv
                      << " selected_count=" << irreps.size() << "\n";
        }
        for (std::size_t k=0; k<irreps.size(); ++k)
        {
            const int assigned_device = (opt.use_multi_gpu && opt.device_count > 1) ? static_cast<int>(k % static_cast<std::size_t>(opt.device_count)) : dev; // v32zu OPT-8: independent irreps can rotate devices without changing output names.
            std::cout << "[L progress] Lbyas=" << par.Lbyas << " irrep " << (k+1) << "/" << irreps.size()
                      << " label=" << irreps[k].label << " device_id=" << assigned_device << "\n";
            irrep_timings.push_back(v32zt_generate_one_irrep(irreps[k], par, E0, E1, coarseN, outdir, opt, assigned_device, debug));
        }

        V32ztTiming sum;
        for (const auto& t : irrep_timings) {
            sum.basis_sec += t.basis_sec;
            sum.f2k2_sec += t.f2k2_sec;
            sum.gf3_build_sec += t.gf3_build_sec;
            sum.solve_sec += t.solve_sec;
            sum.vsel_sec += t.vsel_sec;
            sum.write_sec += t.write_sec;
            sum.rows += t.rows;
            sum.chunks += t.chunks;
            sum.bytes_written += t.bytes_written;
            sum.max_total_dim = std::max(sum.max_total_dim, t.max_total_dim);
            sum.max_vdim = std::max(sum.max_vdim, t.max_vdim);
        }
        const double full_L_sec = Ltimer.elapsed();
        std::cout << "============================================================\n";
        std::cout << "[Lbyas FULL RUN SUMMARY]\n";
        std::cout << "Lbyas                         = " << par.Lbyas << "\n";
        std::cout << "coarseN                       = " << coarseN << "\n";
        std::cout << "default irreps completed       = " << irreps.size() << " / " << (opt.irrep_filter_csv.empty() ? 5 : static_cast<int>(irreps.size())) << "\n";
        std::cout << "total expected grid rows       = " << (static_cast<long long>(irreps.size()) * coarseN) << "\n";
        std::cout << "total rows actually written    = " << sum.rows << "\n";
        std::cout << "total chunks                   = " << sum.chunks << "\n";
        std::cout << "max total dim                  = " << sum.max_total_dim << "\n";
        std::cout << "max vdim                       = " << sum.max_vdim << "\n";
        std::cout << "basis construction time        = " << sum.basis_sec << " sec\n";
        std::cout << "F2/K2 time                     = " << sum.f2k2_sec << " sec\n";
        std::cout << "G/F3 matrix build time         = " << sum.gf3_build_sec << " sec\n";
        std::cout << "matrix inverse/solve time      = " << sum.solve_sec << " sec\n";
        std::cout << "Vsel construction time         = " << sum.vsel_sec << " sec\n";
        std::cout << "binary write time              = " << sum.write_sec << " sec\n";
        std::cout << "full Lbyas wall time           = " << full_L_sec << " sec (" << v32zt_format_seconds(full_L_sec) << ")\n";
        std::cout << "average time per irrep         = " << (full_L_sec / std::max<double>(1.0, double(irreps.size()))) << " sec\n";
        std::cout << "average time per written row   = " << (full_L_sec / std::max<double>(1.0, double(sum.rows))) << " sec\n";
        std::cout << "output directory               = " << outdir << "\n";
        if (opt.validate_fused_f2k2) {                                                   // v32zu OPT-1: final validation summary for smoke logs/report extraction.
            std::cout << "[F2K2 validation] chunks=" << g_v32zu_f2k2_validation.chunks
                      << " max_abs_diff_F2=" << g_v32zu_f2k2_validation.max_abs_diff_F2
                      << " max_rel_diff_F2=" << g_v32zu_f2k2_validation.max_rel_diff_F2
                      << " max_abs_diff_K2inv=" << g_v32zu_f2k2_validation.max_abs_diff_K2inv
                      << " max_rel_diff_K2inv=" << g_v32zu_f2k2_validation.max_rel_diff_K2inv << "\n";
        }
        std::cout << "[DONE] v32zu cache generation completed\n";
        std::cout << "============================================================\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[v32zt ERROR] " << e.what() << "\n";
        return 1;
    }
}
