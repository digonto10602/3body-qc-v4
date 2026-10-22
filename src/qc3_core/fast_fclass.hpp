// fast_fclass.hpp : class moment tables for the F lattice sum (prompts/002 sub-step 7.7, prompts/001 rung 5).
//
// F-block lattice sum for one spectator momentum (fast_point.hpp):  S_ab = sum_{n in cube} Y_a(r) Y_b(r) w(r^2),
//   r = n + nP * prod1(n.nP),  prod1(d) = (d/nP^2) (1/gamma - 1) - xi/gamma,   w = exp(alpha (x^2 - r^2)) / (x^2 - r^2),
// for l <= 1 the harmonics are Y_a = u_a + v_a . r (u = 1, v = 0 for l = 0; u = 0, v = sqrt3 e_i for l = 1), so only the moments
//   T0 = sum w,  T1_i = sum w r_i,  T2_ij = sum w r_i r_j     (w r^2-ball truncated, see below)
// are needed:  S_ab = u_a u_b T0 + u_a v_b.T1 + u_b v_a.T1 + v_a^T T2 v_b.
//
// Exact algebra, no symmetry assumption beyond lattice-cube invariance:
//  (1) r^2 = n^2 + 2 prod1 (n.nP) + prod1^2 nP^2 depends on n only through the integers (n^2, n.nP): a "class" is the set of cube
//      points with equal (n^2, n.nP).  w is constant on a class.  Per class the energy-independent moments are stored:
//      count, mean n, centred second moments C = sum (n - nbar)(n - nbar)^T (from exact integer sums), and e^{-alpha perp2}.
//      Per energy: sum over classes of w * (count, count*rbar, C + count*rbar rbar^T), rbar = nbar + nP prod1 (no cancellation).
//  (2) r^2 = perp2 + rpar^2, perp2 = n^2 - (n.nP)^2/nP^2 (exact rational), rpar = (d + nP^2 prod1)/|nP|.  Classes with equal
//      d = n.nP form a group sorted by perp2; the ellipsoid r^2 <= R^2 is  perp2 <= R^2 - rpar^2  (binary search), the same
//      point set as the direct sum's truncation (so its truncation guarantee is unchanged).  exp(alpha(x^2-r^2)) =
//      exp(alpha (x^2 - rpar^2)) * exp(-alpha perp2): one exp per group, the second factor is in the table.
//  (3) Tables are built for the canonical nP = (c0 >= c1 >= c2 >= 0) of the cubic-group orbit of nP.  A signed permutation R with
//      R nP = c maps the cube to itself and n.nP -> m.c, so the moments of nP are R^T-rotated moments of c.
//  (4) S depends on the spectator only through (c, |p|^2, energy): momenta with equal (c, |p|^2) share one sum (rotated).
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace qc3_fast {

#ifndef QC3_FACC
#define QC3_FACC double
#endif
using Facc = QC3_FACC;

struct FClass {
    double perp2, ep, cnt;
    double nb[3];
    double C[6];   // xx, xy, xz, yy, yz, zz (centred, about nb)
};

struct FGroup { int dot; int begin, end; double dotd; };

struct FTable {
    int c[3] = {0, 0, 0};
    int nPsq = 0;
    int nmax = 0;
    double alpha = 0.0;
    std::vector<FClass> cl;
    std::vector<FGroup> grp;         // ascending dot
    long long total_points = 0;
};

inline std::shared_ptr<const FTable> build_ftable(const int c[3], int nmax, double alpha) {
    auto t = std::make_shared<FTable>();
    t->c[0] = c[0]; t->c[1] = c[1]; t->c[2] = c[2];
    t->nPsq = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    t->nmax = nmax; t->alpha = alpha;
    const int w = 2 * nmax + 1;
    const int n2max = 3 * nmax * nmax;
    const int dmax = nmax * (c[0] + c[1] + c[2]);
    const long long dense = (long long)(2 * dmax + 1) * (n2max + 1);
    struct Acc { long long cnt = 0, s[3] = {0, 0, 0}, q[6] = {0, 0, 0, 0, 0, 0}; int n2 = 0, dot = 0; };
    std::vector<Acc> acc;
    std::unordered_map<long long, int> sparse;
    std::vector<int> dense_idx;
    const bool use_dense = dense <= 6000000;
    if (use_dense) dense_idx.assign(std::size_t(dense), -1);
    else sparse.reserve(40000);
    for (int ix = -nmax; ix <= nmax; ++ix)
        for (int iy = -nmax; iy <= nmax; ++iy)
            for (int iz = -nmax; iz <= nmax; ++iz) {
                const int n2 = ix * ix + iy * iy + iz * iz;
                const int dot = ix * c[0] + iy * c[1] + iz * c[2];
                int* slot;
                if (use_dense) slot = &dense_idx[std::size_t((long long)(dot + dmax) * (n2max + 1) + n2)];
                else {
                    auto it = sparse.emplace((long long)(dot + dmax) * (n2max + 1) + n2, -1).first;
                    slot = &it->second;
                }
                if (*slot < 0) { *slot = int(acc.size()); acc.emplace_back(); acc.back().n2 = n2; acc.back().dot = dot; }
                Acc& a = acc[std::size_t(*slot)];
                a.cnt += 1;
                a.s[0] += ix; a.s[1] += iy; a.s[2] += iz;
                a.q[0] += ix * ix; a.q[1] += ix * iy; a.q[2] += ix * iz; a.q[3] += iy * iy; a.q[4] += iy * iz; a.q[5] += iz * iz;
            }
    (void)w;
    std::vector<int> order(acc.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (acc[std::size_t(a)].dot != acc[std::size_t(b)].dot) return acc[std::size_t(a)].dot < acc[std::size_t(b)].dot;
        return acc[std::size_t(a)].n2 < acc[std::size_t(b)].n2;
    });
    t->cl.reserve(acc.size());
    for (int idx : order) {
        const Acc& a = acc[std::size_t(idx)];
        FClass k;
        const long long numer = (t->nPsq > 0) ? ((long long)a.n2 * t->nPsq - (long long)a.dot * a.dot) : 0;
        k.perp2 = (t->nPsq > 0) ? double(numer) / double(t->nPsq) : double(a.n2);
        k.ep = std::exp(-alpha * k.perp2);
        k.cnt = double(a.cnt);
        for (int i = 0; i < 3; ++i) k.nb[i] = double(a.s[i]) / double(a.cnt);
        const int ii[6] = {0, 0, 0, 1, 1, 2}, jj[6] = {0, 1, 2, 1, 2, 2};
        for (int m = 0; m < 6; ++m)   // exact integer numerator of the centred second moment
            k.C[m] = double(a.cnt * a.q[m] - a.s[ii[m]] * a.s[jj[m]]) / double(a.cnt);
        t->total_points += a.cnt;
        if (t->grp.empty() || t->grp.back().dot != a.dot) {
            FGroup g; g.dot = a.dot; g.begin = int(t->cl.size()); g.end = g.begin; g.dotd = double(a.dot);
            t->grp.push_back(g);
        }
        t->cl.push_back(k);
        t->grp.back().end = int(t->cl.size());
    }
    return t;
}

struct FTableCache {
    std::mutex mtx;
    std::unordered_map<std::string, std::shared_ptr<const FTable>> map;
    long long hits = 0, misses = 0;
};
inline FTableCache& ftable_cache() { static FTableCache c; return c; }

inline void clear_ftable_cache() {
    FTableCache& c = ftable_cache();
    std::lock_guard<std::mutex> lock(c.mtx);
    c.map.clear();
}
inline long long ftable_entries() {
    FTableCache& c = ftable_cache();
    std::lock_guard<std::mutex> lock(c.mtx);
    return (long long)c.map.size();
}

inline std::shared_ptr<const FTable> ftable_cached(const int c[3], int nmax, double alpha) {
    FTableCache& cache = ftable_cache();
    char key[96];
    std::snprintf(key, sizeof key, "%d,%d,%d,%d,%.17g", c[0], c[1], c[2], nmax, alpha);
    {
        std::lock_guard<std::mutex> lock(cache.mtx);
        auto it = cache.map.find(key);
        if (it != cache.map.end()) { ++cache.hits; return it->second; }
        ++cache.misses;
    }
    auto t = build_ftable(c, nmax, alpha);
    std::lock_guard<std::mutex> lock(cache.mtx);
    cache.map.emplace(key, t);
    return t;
}

// Moments of one canonical (c, |p|^2) sum.
struct FMoments {
    double T0 = 0.0, T1[3] = {0, 0, 0}, T2[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    long long terms = 0, classes = 0;
};

// Moments of the ellipsoid-truncated sum for the canonical nP = c:  w = exp(alpha (x2 - r^2)) / (x2 - r^2),  r^2 <= R2.
inline FMoments class_moments(const FTable& t, double x2, double gamma, double xii, double R2) {
    const double alpha = t.alpha;
    const bool boosted = (t.nPsq > 0);
    const double nPsq = double(t.nPsq), nPn = std::sqrt(nPsq);
    const double xibygamma = xii / gamma;
    const double invg_m1 = 1.0 / gamma - 1.0;
    const double cv[3] = {double(t.c[0]), double(t.c[1]), double(t.c[2])};
    Facc T0 = 0, T1[3] = {0, 0, 0}, T2[6] = {0, 0, 0, 0, 0, 0};
    long long terms = 0, classes = 0;
    for (const FGroup& g : t.grp) {
        double prod1 = 0.0, rpar2 = 0.0;
        if (boosted) {
            prod1 = (g.dotd / nPsq) * invg_m1 - xibygamma;
            const double rpar = (g.dotd + nPsq * prod1) / nPn;
            rpar2 = rpar * rpar;
        }
        const double avail = R2 - rpar2;
        if (avail < 0.0) continue;
        const FClass* b = t.cl.data() + g.begin;
        const FClass* e = t.cl.data() + g.end;
        // classes are sorted by perp2 ascending: keep perp2 <= avail
        const FClass* hi = std::upper_bound(b, e, avail, [](double v, const FClass& k) { return v < k.perp2; });
        if (hi == b) continue;
        const double G = x2 - rpar2;
        const double off[3] = {cv[0] * prod1, cv[1] * prod1, cv[2] * prod1};
        const double eg = std::exp(alpha * G);
        // ep = exp(-alpha perp2) >= exp(-alpha avail) for every class used here: no underflow if alpha*avail < 575;
        // eg*ep is the exp of the direct sum.  Otherwise (huge R^2 or alpha) fall back to one exp per class.
        const bool split = (alpha * avail < 575.0) && std::isfinite(eg) && eg > 0.0 && eg < 1e300;
        for (const FClass* k = b; k != hi; ++k) {
            const double dd = G - k->perp2;
            const double w = (split ? eg * k->ep : std::exp(alpha * dd)) / dd;
            const double wc = w * k->cnt;
            const double rx = k->nb[0] + off[0], ry = k->nb[1] + off[1], rz = k->nb[2] + off[2];
            const double wx = wc * rx, wy = wc * ry, wz = wc * rz;
            T0 += wc;
            T1[0] += wx; T1[1] += wy; T1[2] += wz;
            T2[0] += wx * rx + w * k->C[0];
            T2[1] += wx * ry + w * k->C[1];
            T2[2] += wx * rz + w * k->C[2];
            T2[3] += wy * ry + w * k->C[3];
            T2[4] += wy * rz + w * k->C[4];
            T2[5] += wz * rz + w * k->C[5];
            terms += (long long)k->cnt;
        }
        classes += hi - b;
    }
    FMoments m;
    m.T0 = double(T0);
    for (int i = 0; i < 3; ++i) m.T1[i] = double(T1[i]);
    m.T2[0][0] = double(T2[0]); m.T2[0][1] = m.T2[1][0] = double(T2[1]); m.T2[0][2] = m.T2[2][0] = double(T2[2]);
    m.T2[1][1] = double(T2[3]); m.T2[1][2] = m.T2[2][1] = double(T2[4]); m.T2[2][2] = double(T2[5]);
    m.terms = terms; m.classes = classes;
    return m;
}

}  // namespace qc3_fast
