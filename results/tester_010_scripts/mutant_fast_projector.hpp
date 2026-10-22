// fast_projector.hpp : the fast path's OWN in-memory projector cache (prompts/002 section 2.6, sub-step 7.0).
//
// The legacy cache (g_projector_cache_v30q, keyed by basis_signature_v30q) omits eig_tol, norm_tol, proj_tol (report 006
// note 2: eig_tol=1e-30 after a default call returned the 17 default columns instead of 5).  The legacy key is NOT edited
// (oracle include chain; the oracle clears that map on every call).  Here the key is
//     basis_signature_v30q (irrep, parity, nnP, full (n,l,m) lists of both channels)
//   + waves_1, waves_2 + eig_tol, norm_tol, proj_tol (17 significant digits) + projector method.
// Method "eigen": the EXISTING build_cached_projector_v30q builder, cached here (cross-check method).
// Method "orbit" (prompts/002 section 2.1, sub-step 7.2): V_sel from the orbits of the little group, no dense P, no
//   eigen-solve.  Conventions are those of the oracle (convention c0 of fv_projector_cartesian_l1_v30m.hpp):
//     U(R) e_(k,l,m) = Pi(R) sum_m' D^(l)_{m'm}(R) e_(Rk,l,m'),   P = d_Gamma/|LG| sum_R chi_Gamma(R) U(R),
//   using the SAME little_group / chi / cubic_transf / parity_factor_for_l / D_real_element_v30m functions.
//   For a 1-dim irrep, P U(R') = chi(R') P, so Im P restricted to a (flavour, orbit, l) block is spanned by
//   P e_(k0,l,m), m = -l..l, k0 = first momentum of the orbit in basis order.  Those <= 2l+1 sparse vectors are
//   orthonormalised by an SVD (left singular vectors with s > 1e-8 s_max).  For d_Gamma > 1 (not in production,
//   implemented but covered only by the self-check) the candidates are all rows of the block.
//   Mandatory self-checks (throw "fast: PROJECTOR_SELF_CHECK"): ||P V - V||_F, ||V^T V - 1||_F and |tr P - Nproj|
//   (tr P = rank P for a projector, so the last one detects an incomplete candidate set).
#pragma once
#include <random>
#include <cstdlib>


#include <mutex>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>

#include "oracle_point.hpp"

namespace qc3_fast {

struct FastProjector {
    qc3_oracle::Mat Vsel;      // N x Nproj (complex, as the legacy builder returns it)
    int Nproj = 0;
    std::string method;
};

struct ProjectorCache {
    std::mutex mtx;
    std::unordered_map<std::string, FastProjector> map;
    long long hits = 0, misses = 0;
};

inline ProjectorCache& projector_cache() {
    static ProjectorCache c;
    return c;
}

inline void clear_fast_caches() {
    ProjectorCache& c = projector_cache();
    std::lock_guard<std::mutex> lock(c.mtx);
    c.map.clear();
}

// ---------------------------------------------------------------------------------------------------------------
// Orbit/character construction (method "orbit")
// ---------------------------------------------------------------------------------------------------------------
namespace orbit_detail {

inline long long row_key(int fl, int nx, int ny, int nz, int ell, int m) {
    if (std::abs(nx) > 63 || std::abs(ny) > 63 || std::abs(nz) > 63 || ell < 0 || ell > 15 || std::abs(m) > 15 || fl < 0 || fl > 3)
        throw std::runtime_error("fast: PROJECTOR basis row outside the key range");
    long long k = fl;
    k = k * 128 + (nx + 64); k = k * 128 + (ny + 64); k = k * 128 + (nz + 64);
    k = k * 16 + ell; k = k * 32 + (m + 16);
    return k;
}

struct Entry { int tgt; double w; };   // U(R) e_src has component w (= chi Pi D, an integer) at row tgt

}  // namespace orbit_detail

inline Eigen::MatrixXd build_orbit_V(const qc3_oracle::Kin& k) {
    using namespace orbit_detail;
    const auto& basis = k.basis;
    const int N = int(basis.size());
    std::unordered_map<long long, int> index;
    index.reserve(std::size_t(2 * N));
    for (int j = 0; j < N; ++j) {
        const auto& b = basis[std::size_t(j)];
        index.emplace(row_key(b[0], b[1], b[2], b[3], b[4], b[5]), j);
    }
    if (int(index.size()) != N) throw std::runtime_error("fast: PROJECTOR duplicate basis rows");

    const std::string& irrep = k.spec.irrep;
    const std::vector<std::vector<int>> LG = little_group(k.nnP_config);
    const int dI = fvproj_v30j::irrep_dimension(irrep);
    const double scale = double(dI) / double(LG.size());

    // sparse U(R) weighted by chi(R): columns grouped by source row
    std::vector<std::vector<Entry>> bysrc(static_cast<std::size_t>(N));
    long long trace_int = 0;   // sum_R chi(R) tr U(R)  (an integer)
    for (const auto& R : LG) {
        const double chiR = double(chi(R, irrep, k.nnP_config));
        const double pf = fvproj_v30j::parity_factor_for_l(k.par.parity, R, 0, 1);   // parity_mode 1: Pi(R), l-independent
        for (int j = 0; j < N; ++j) {
            const auto& b = basis[std::size_t(j)];
            const std::vector<int> Rn = cubic_transf({b[1], b[2], b[3]}, R);
            const int ell = b[4], ms = b[5];
            for (int mt = -ell; mt <= ell; ++mt) {
                const double d = fvproj_v30j::D_real_element_v30m(ell, mt, ms, R);
                if (d == 0.0) continue;
                auto it = index.find(row_key(b[0], Rn[0], Rn[1], Rn[2], ell, mt));
                if (it == index.end())
                    throw std::runtime_error("fast: PROJECTOR_SELF_CHECK basis is not closed under the little group");
                bysrc[std::size_t(j)].push_back({it->second, chiR * pf * d});
                if (it->second == j) trace_int += llround(chiR * pf * d);
            }
        }
    }
    const double trP = scale * double(trace_int);

    // blocks: (flavour, orbit, l) in order of first appearance in the basis
    std::vector<std::vector<double>> cols;   // each: dense length-N column (orthonormal after the SVD)
    std::unordered_set<long long> seen;
    for (int j0 = 0; j0 < N; ++j0) {
        const auto& b0 = basis[std::size_t(j0)];
        const int fl = b0[0], ell = b0[4];
        if (seen.count(row_key(fl, b0[1], b0[2], b0[3], ell, 0))) continue;
        std::vector<int> cand;   // candidate source rows
        std::vector<std::vector<int>> orbit;
        for (const auto& R : LG) {
            std::vector<int> Rn = cubic_transf({b0[1], b0[2], b0[3]}, R);
            bool dup = false;
            for (const auto& o : orbit) if (o == Rn) { dup = true; break; }
            if (!dup) orbit.push_back(Rn);
        }
        for (const auto& o : orbit) seen.insert(row_key(fl, o[0], o[1], o[2], ell, 0));
        for (const auto& o : orbit) {
            if (dI == 1 && o != std::vector<int>{b0[1], b0[2], b0[3]}) continue;   // 1-dim irrep: representative only
            for (int m = -ell; m <= ell; ++m) {
                auto it = index.find(row_key(fl, o[0], o[1], o[2], ell, m));
                if (it != index.end()) cand.push_back(it->second);
            }
        }
        // P e_j for the candidates (integer accumulation, scaled afterwards), on the union of their supports
        std::vector<int> support;
        for (int j : cand) for (const auto& e : bysrc[std::size_t(j)]) support.push_back(e.tgt);
        std::sort(support.begin(), support.end());
        support.erase(std::unique(support.begin(), support.end()), support.end());
        if (support.empty() || cand.empty()) continue;
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(long(support.size()), long(cand.size()));
        for (std::size_t c = 0; c < cand.size(); ++c)
            for (const auto& e : bysrc[std::size_t(cand[c])]) {
                const long r = std::lower_bound(support.begin(), support.end(), e.tgt) - support.begin();
                C(r, long(c)) += e.w;
            }
        C *= scale;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeThinU);
        const Eigen::VectorXd& s = svd.singularValues();
        const double smax = s.size() ? s(0) : 0.0;
        if (!(smax > 1e-8)) continue;   // block has no component in the irrep (entries are exact multiples of 1/|LG|)
        for (long i = 0; i < s.size(); ++i) {
            if (!(s(i) > 1e-8 * smax)) break;
            std::vector<double> col(std::size_t(N), 0.0);
            for (std::size_t r = 0; r < support.size(); ++r) col[std::size_t(support[r])] = svd.matrixU()(long(r), i);
            cols.push_back(std::move(col));
        }
    }
    const int Np = int(cols.size());
    Eigen::MatrixXd V(N, Np);
    for (int c = 0; c < Np; ++c) for (int r = 0; r < N; ++r) V(r, c) = cols[std::size_t(c)][std::size_t(r)];

    // ---- mandatory self-checks ----
    if (Np > 0) {
        Eigen::MatrixXd PV = Eigen::MatrixXd::Zero(N, Np);
        for (int j = 0; j < N; ++j)
            for (const auto& e : bysrc[std::size_t(j)]) PV.row(e.tgt) += e.w * V.row(j);
        PV *= scale;
        const double res = (PV - V).norm();
        const double orth = (V.transpose() * V - Eigen::MatrixXd::Identity(Np, Np)).norm();
        if (!(res <= 1e-10 * std::sqrt(double(Np))))
            throw std::runtime_error("fast: PROJECTOR_SELF_CHECK ||P V - V|| = " + std::to_string(res));
        if (!(orth <= 1e-12))
            throw std::runtime_error("fast: PROJECTOR_SELF_CHECK ||V^T V - 1|| = " + std::to_string(orth));
    }
    if (!(std::abs(trP - double(Np)) <= 1e-9))
        throw std::runtime_error("fast: PROJECTOR_SELF_CHECK tr P = " + std::to_string(trP) + " but Nproj = " + std::to_string(Np));
    return V;
}

inline std::string projector_key(const qc3_oracle::Kin& k, const std::string& method) {
    std::ostringstream os;
    os << std::setprecision(17);
    os << basis_signature_v30q(k.plm_config, k.np_config, k.klm_config, k.nk_config, k.spec.irrep, k.nnP_config, k.par.parity);
    os << ";w1=";
    for (int w : k.par.waves_vec_1) os << w << ',';
    os << ";w2=";
    for (int w : k.par.waves_vec_2) os << w << ',';
    os << ";eig_tol=" << k.par.eig_tol << ";norm_tol=" << k.par.norm_tol << ";proj_tol=" << k.par.proj_tol
       << ";method=" << method;
    return os.str();
}

// Returns a copy of the cached projector (building it on a miss). `hit` tells whether the entry existed.
inline FastProjector fast_projector_cached(const qc3_oracle::Kin& k, const std::string& method, bool& hit) {
    ProjectorCache& c = projector_cache();
    const std::string key = projector_key(k, method);
    {
        std::lock_guard<std::mutex> lock(c.mtx);
        auto it = c.map.find(key);
        if (it != c.map.end()) { ++c.hits; hit = true; return it->second; }
    }
    hit = false;
    FastProjector fp;
    fp.method = method;
    if (method == "eigen") {
        CachedProjectorV30q built = build_cached_projector_v30q(k.plm_config, k.np_config, k.klm_config, k.nk_config,
                                                                k.spec.irrep, k.nnP_config, k.par);
        fp.Vsel = built.Vsel;
        fp.Nproj = int(built.Vsel.cols());
    } else if (method == "orbit") {
        fp.Vsel = build_orbit_V(k).cast<comp>();
        fp.Nproj = int(fp.Vsel.cols());
        if (std::getenv("TESTER_ROT")) {   // TESTER MUTANT: V -> V U, U a fixed random complex unitary (exercises the complex-V tail branch)
            const int n = fp.Nproj;
            std::mt19937_64 g(20261022); std::normal_distribution<double> nd(0.0, 1.0);
            Eigen::MatrixXcd A(n, n); for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) A(i, j) = comp(nd(g), nd(g));
            Eigen::HouseholderQR<Eigen::MatrixXcd> qr(A); Eigen::MatrixXcd U = qr.householderQ() * Eigen::MatrixXcd::Identity(n, n);
            fp.Vsel = (fp.Vsel * U).eval();
        }
    } else {
        throw std::runtime_error("fast projector method '" + method + "' is unknown (expected 'orbit' or 'eigen')");
    }
    std::lock_guard<std::mutex> lock(c.mtx);
    ++c.misses;
    c.map.emplace(key, fp);
    return fp;
}

}  // namespace qc3_fast
