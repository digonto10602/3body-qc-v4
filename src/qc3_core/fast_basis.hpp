// fast_basis.hpp : threshold basis builder and dimension jumps (prompts/002 section 2.2, sub-step 7.3).
//
// The oracle's basis builder config_maker_4_momentum_first scans all (2*nmax+1)^3 = 41^3 candidate spectator
// momenta n for EVERY energy and keeps n iff
//     En.real() > En_min_plus_for_config(p, P, mi, mj, mk, 0.02, epsilon_h).real(),   p = (2 pi / L) n.
// The right-hand side does not depend on the energy.  Here it is computed once per (sector, channel) with the
// SAME function and the SAME arguments (same comp arithmetic) and stored in loop order (nx outer ... nz inner);
// per energy, n is kept iff E_lab > t[n] -- the identical double comparison, so the selection and the order are
// bit-identical to the oracle's by construction.  The (ell, m) list is expanded exactly as the oracle does.
// The sorted distinct thresholds, converted to E* = sqrt(t^2 - P^2), are the dimension jumps of the basis.
#pragma once

#include <algorithm>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <unordered_map>

#include "oracle_point.hpp"

namespace qc3_fast {

// Parameters, spectrum label and total momentum: the first half of qc3_oracle::make_kinematics (identical
// statements; the oracle file is not edited), without the basis.
inline qc3_oracle::Kin make_setup(const qc3_oracle::Input& in) {
    using namespace k3df_fit_v32f;
    qc3_oracle::Kin k;
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
    k.par = make_base_physics(s);
    k.par.waves_vec_1 = in.waves_1;
    k.par.waves_vec_2 = in.waves_2;
    k.spec = parse_label(in.irrep);

    const comp pi = std::acos(-1.0);
    k.L = k.par.L();
    const comp twopibyL = comp(2.0, 0.0) * pi / comp(k.L, 0.0);
    k.total_P.assign(3, comp(0.0, 0.0)); k.nnP_config.assign(3, comp(0.0, 0.0));
    for (int a = 0; a < 3; ++a) {
        k.total_P[a] = twopibyL * double(k.spec.nnP[a]);
        k.nnP_config[a] = comp(k.spec.nnP[a], 0.0);
    }
    return k;
}

// ---- threshold tables: one per (sector, channel), in memory ----
struct ThresholdCache {
    std::mutex mtx;
    std::unordered_map<std::string, std::vector<double>> map;   // node-based: references stay valid
};

inline ThresholdCache& threshold_cache() {
    static ThresholdCache c;
    return c;
}

inline void clear_threshold_tables() {
    ThresholdCache& c = threshold_cache();
    std::lock_guard<std::mutex> lock(c.mtx);
    c.map.clear();
}

inline long long threshold_entries() {
    ThresholdCache& c = threshold_cache();
    std::lock_guard<std::mutex> lock(c.mtx);
    return (long long)c.map.size();
}

// t[n] = En_min_plus_for_config(p, P, mi, mj, mk, 0.02, epsilon_h).real() in the loop order of the oracle.
inline std::vector<double> compute_thresholds(const qc3_oracle::Kin& k, double mi, double mj, double mk) {
    const comp pi = std::acos(-1.0);
    const comp twopibyL = ((comp)2.0 * pi) / ((comp)k.L);
    const int nmax = static_cast<int>(k.par.max_shell_num);
    const int side = 2 * nmax + 1;
    std::vector<double> t;
    t.reserve(std::size_t(side) * std::size_t(side) * std::size_t(side));
    for (int i = -nmax; i < nmax + 1; ++i)
        for (int j = -nmax; j < nmax + 1; ++j)
            for (int l = -nmax; l < nmax + 1; ++l) {
                std::vector<comp> pvec(3);
                pvec[0] = twopibyL * ((comp)i);
                pvec[1] = twopibyL * ((comp)j);
                pvec[2] = twopibyL * ((comp)l);
                t.push_back(En_min_plus_for_config(pvec, k.total_P, mi, mj, mk, 0.02, k.par.epsilon_h).real());
            }
    return t;
}

// channel 1: spectator K (mi = mK, pair K pi); channel 2: spectator pi (mi = mpi, pair K K) -- as make_kinematics.
inline const std::vector<double>& threshold_table(const qc3_oracle::Kin& k, int channel) {
    const double mi = (channel == 1) ? k.par.atmK : k.par.atmpi;
    const double mj = k.par.atmK;
    const double mk = (channel == 1) ? k.par.atmpi : k.par.atmK;
    std::ostringstream os;
    os << std::setprecision(17) << "P=" << k.spec.nnP[0] << ',' << k.spec.nnP[1] << ',' << k.spec.nnP[2] << ";L=" << k.L
       << ";eps_h=" << k.par.epsilon_h << ";nmax=" << k.par.max_shell_num << ";m=" << mi << ',' << mj << ',' << mk;
    const std::string key = os.str();
    ThresholdCache& c = threshold_cache();
    {
        std::lock_guard<std::mutex> lock(c.mtx);
        auto it = c.map.find(key);
        if (it != c.map.end()) return it->second;
    }
    std::vector<double> t = compute_thresholds(k, mi, mj, mk);
    std::lock_guard<std::mutex> lock(c.mtx);
    return c.map.emplace(key, std::move(t)).first->second;
}

// One channel's configuration, bit-identical to config_maker_4_momentum_first (plm: comp momenta + ell, m; n: ints).
inline void fast_config(const std::vector<double>& thr, double E_lab, const std::vector<int>& waves, const qc3_oracle::Kin& k,
                        std::vector<std::vector<comp>>& plm_config, std::vector<std::vector<int>>& n_config) {
    const comp pi = std::acos(-1.0);
    const comp twopibyL = ((comp)2.0 * pi) / ((comp)k.L);
    const int nmax = static_cast<int>(k.par.max_shell_num);
    std::vector<std::vector<int>> ell_m_vec(2, std::vector<int>());
    for (int iw = 0; iw < static_cast<int>(waves.size()); ++iw) {
        const int ell = waves[std::size_t(iw)];
        std::vector<int> m_vec;
        ell_m_vector(ell, m_vec);
        for (int jm = 0; jm < static_cast<int>(m_vec.size()); ++jm) {
            ell_m_vec[0].push_back(ell);
            ell_m_vec[1].push_back(m_vec[std::size_t(jm)]);
        }
    }
    std::size_t nsel = 0;
    for (double t : thr) if (E_lab > t) ++nsel;
    const std::size_t rows = nsel * ell_m_vec[0].size();
    for (int a = 0; a < 5; ++a) { plm_config[std::size_t(a)].reserve(rows); n_config[std::size_t(a)].reserve(rows); }

    std::size_t idx = 0;
    for (int i = -nmax; i < nmax + 1; ++i)
        for (int j = -nmax; j < nmax + 1; ++j)
            for (int l = -nmax; l < nmax + 1; ++l, ++idx) {
                if (!(E_lab > thr[idx])) continue;
                const comp px = twopibyL * ((comp)i);
                const comp py = twopibyL * ((comp)j);
                const comp pz = twopibyL * ((comp)l);
                for (int lm = 0; lm < static_cast<int>(ell_m_vec[0].size()); ++lm) {
                    const int ell = ell_m_vec[0][std::size_t(lm)];
                    const int proj_m = ell_m_vec[1][std::size_t(lm)];
                    plm_config[0].push_back(px);
                    plm_config[1].push_back(py);
                    plm_config[2].push_back(pz);
                    plm_config[3].push_back(ell);
                    plm_config[4].push_back(proj_m);
                    n_config[0].push_back(i);
                    n_config[1].push_back(j);
                    n_config[2].push_back(l);
                    n_config[3].push_back(ell);
                    n_config[4].push_back(proj_m);
                }
            }
}

// Same Kin as qc3_oracle::make_kinematics (the oracle keeps config_maker_4_momentum_first).
inline qc3_oracle::Kin fast_kinematics(const qc3_oracle::Input& in) {
    qc3_oracle::Kin k = make_setup(in);
    const comp Ecm_c(in.Ecm, 0.0);
    k.En_c = Ecm_to_E(Ecm_c, k.total_P);
    k.plm_config.assign(5, {}); k.klm_config.assign(5, {});
    k.np_config.assign(5, {}); k.nk_config.assign(5, {});
    const double E_lab = k.En_c.real();
    fast_config(threshold_table(k, 1), E_lab, k.par.waves_vec_1, k, k.plm_config, k.np_config);
    fast_config(threshold_table(k, 2), E_lab, k.par.waves_vec_2, k, k.klm_config, k.nk_config);
    k.N1 = int(k.plm_config[0].size());
    k.N2 = int(k.klm_config[0].size());
    k.N = k.N1 + k.N2;
    if (k.N <= 0) throw std::runtime_error("EMPTY_CONFIG (no basis states at this energy)");
    k.basis.reserve(std::size_t(k.N));
    for (int a = 0; a < k.N1; ++a)
        k.basis.push_back({1, k.np_config[0][a], k.np_config[1][a], k.np_config[2][a], k.np_config[3][a], k.np_config[4][a]});
    for (int a = 0; a < k.N2; ++a)
        k.basis.push_back({2, k.nk_config[0][a], k.nk_config[1][a], k.nk_config[2][a], k.nk_config[3][a], k.nk_config[4][a]});
    return k;
}

// Sorted E_cm strictly inside (E_lo, E_hi) at which the basis changes (either channel).  Thresholds are E_lab values;
// E* = sqrt(t^2 - P^2).  Thresholds equal within 1e-12 relative (orbit partners: round-off only) are one jump.
inline std::vector<double> dimension_jumps(const qc3_oracle::Input& in, double E_lo, double E_hi) {
    const qc3_oracle::Kin k = make_setup(in);
    double Psq = 0.0;
    for (int a = 0; a < 3; ++a) Psq += k.total_P[std::size_t(a)].real() * k.total_P[std::size_t(a)].real();
    std::vector<double> inside;
    for (int ch = 1; ch <= 2; ++ch) {
        for (double t : threshold_table(k, ch)) {
            const double Es = std::sqrt(t * t - Psq);
            if (Es > E_lo && Es < E_hi) inside.push_back(t);
        }
    }
    std::sort(inside.begin(), inside.end());
    std::vector<double> out;
    double last = 0.0;
    bool have = false;
    for (double t : inside) {
        if (have && std::abs(t - last) <= 1e-12 * t) continue;
        out.push_back(std::sqrt(t * t - Psq));
        last = t; have = true;
    }
    return out;
}

}  // namespace qc3_fast
