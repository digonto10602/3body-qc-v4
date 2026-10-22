// tester harness (not production): slogdet_lu of fast_point.hpp on known matrices; singular / non-finite must throw.
#include "fast_point.hpp"
#include <iostream>
#include <random>
using namespace qc3_fast;
int main() {
    std::mt19937_64 g(20261022);
    std::normal_distribution<double> nd(0.0, 1.0);
    // 1. random real symmetric matrices of several sizes: compare with log|det| from eigenvalues (SelfAdjointEigenSolver) and sign
    double worst_ln = 0; int sign_bad = 0; int nneg = 0;
    for (int n : {1, 2, 5, 17, 40, 103}) for (int rep = 0; rep < 30; ++rep) {
        Eigen::MatrixXd A(n, n);
        for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j) { A(i, j) = nd(g); A(j, i) = A(i, j); }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
        double ln = 0; int neg = 0;
        for (int k = 0; k < n; ++k) { ln += std::log(std::abs(es.eigenvalues()[k])); if (es.eigenvalues()[k] < 0) ++neg; }
        int sgn = (neg % 2 == 0) ? 1 : -1; if (sgn < 0) ++nneg;
        int s; double si, l; slogdet_lu(A, s, si, l);
        worst_ln = std::max(worst_ln, std::abs(l - ln) / std::max(1.0, std::abs(ln)));
        if (s != sgn) ++sign_bad;
    }
    std::cout << "real: worst rel ln|det| diff " << worst_ln << ", sign mismatches " << sign_bad << " (negative-det cases: " << nneg << ")\n";
    // 2. complex Hermitian
    worst_ln = 0; sign_bad = 0; double worst_si = 0;
    for (int n : {2, 7, 33}) for (int rep = 0; rep < 20; ++rep) {
        Eigen::MatrixXcd A(n, n);
        for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j) { comp z(nd(g), i == j ? 0.0 : nd(g)); A(i, j) = z; A(j, i) = std::conj(z); }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
        double ln = 0; int neg = 0;
        for (int k = 0; k < n; ++k) { ln += std::log(std::abs(es.eigenvalues()[k])); if (es.eigenvalues()[k] < 0) ++neg; }
        int sgn = (neg % 2 == 0) ? 1 : -1;
        int s; double si, l; slogdet_lu(A, s, si, l);
        worst_ln = std::max(worst_ln, std::abs(l - ln) / std::max(1.0, std::abs(ln))); worst_si = std::max(worst_si, si);
        if (s != sgn) ++sign_bad;
    }
    std::cout << "complex Hermitian: worst rel ln|det| diff " << worst_ln << ", sign mismatches " << sign_bad << ", max |Im phase| " << worst_si << "\n";
    // 3. exact singular (zero row) -> throw ; non-finite -> throw
    int s; double si, l;
    { Eigen::MatrixXd A = Eigen::MatrixXd::Identity(4, 4); A.row(2).setZero();
      try { slogdet_lu(A, s, si, l); std::cout << "singular: NO THROW (l=" << l << ")\n"; } catch (const std::exception& e) { std::cout << "singular: THROW " << e.what() << "\n"; } }
    { Eigen::MatrixXd A = Eigen::MatrixXd::Identity(4, 4); A(1, 1) = std::nan("");
      try { slogdet_lu(A, s, si, l); std::cout << "nan: NO THROW\n"; } catch (const std::exception& e) { std::cout << "nan: THROW " << e.what() << "\n"; } }
    { Eigen::MatrixXd A = Eigen::MatrixXd::Identity(4, 4); A(3, 0) = INFINITY;
      try { slogdet_lu(A, s, si, l); std::cout << "inf: NO THROW\n"; } catch (const std::exception& e) { std::cout << "inf: THROW " << e.what() << "\n"; } }
    { Eigen::MatrixXd A(2, 2); A << 0, 1, 1, 0;   // needs a row swap: det = -1
      slogdet_lu(A, s, si, l); std::cout << "swap matrix det=-1: sign " << s << " ln " << l << "\n"; }
    { Eigen::MatrixXd A = Eigen::MatrixXd::Identity(3, 3) * 1e-200;   // tiny: ln|det| = 3 ln(1e-200)
      slogdet_lu(A, s, si, l); std::cout << "tiny 1e-200 I3: sign " << s << " ln " << l << " expected " << 3 * std::log(1e-200) << "\n"; }
    return 0;
}
