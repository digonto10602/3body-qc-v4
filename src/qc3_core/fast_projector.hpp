// fast_projector.hpp : the fast path's OWN in-memory projector cache (prompts/002 section 2.6, sub-step 7.0).
//
// The legacy cache (g_projector_cache_v30q, keyed by basis_signature_v30q) omits eig_tol, norm_tol, proj_tol (report 006
// note 2: eig_tol=1e-30 after a default call returned the 17 default columns instead of 5).  The legacy key is NOT edited
// (oracle include chain; the oracle clears that map on every call).  Here the key is
//     basis_signature_v30q (irrep, parity, nnP, full (n,l,m) lists of both channels)
//   + waves_1, waves_2 + eig_tol, norm_tol, proj_tol (17 significant digits) + projector method.
// Method "eigen": the EXISTING build_cached_projector_v30q builder, cached here.  ("orbit": prompts/002 step 7.2.)
#pragma once

#include <mutex>
#include <sstream>
#include <iomanip>
#include <unordered_map>

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
    } else {
        throw std::runtime_error("fast projector method '" + method + "' is not implemented (prompts/002 step 7.2)");
    }
    std::lock_guard<std::mutex> lock(c.mtx);
    ++c.misses;
    c.map.emplace(key, fp);
    return fp;
}

}  // namespace qc3_fast
