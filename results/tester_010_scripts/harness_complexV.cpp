// tester harness: complex-V branch of qc3_fast::evaluate.  The mutant fast_projector.hpp (shadowing the production one by -I order) rotates V -> V U.
#include "fast_point.hpp"
#include <iostream>
#include <random>
using namespace qc3_fast;
int main() {
    for (auto pr : std::vector<std::pair<std::string, double>>{{"100_A2", 0.2915561}, {"111_A2", 0.32036}, {"110_A2", 0.30733}, {"000_A1m", 0.34538}, {"200_A2", 0.34738}}) {
        qc3_oracle::Input in; in.irrep = pr.first; in.Ecm = pr.second; in.k2_pwave_form = "coded_v2";
        in.c = {73735.840894011912, -972421.14060757787, 347174.05548116949, -1226756.7068845264};
        qc3_fast::Options opt; opt.want_F3inv = false; opt.want_F = opt.want_G = opt.want_K2inv = opt.want_F3 = opt.want_Vsel = false;
        clear_fast_caches();
        unsetenv("TESTER_ROT");
        auto o1 = qc3_fast::evaluate(in, opt);
        clear_fast_caches();
        setenv("TESTER_ROT", "1", 1);
        auto o2 = qc3_fast::evaluate(in, opt);
        // reconstruct U exactly as the mutant does
        const int n = o1.Nproj;
        std::mt19937_64 g(20261022); std::normal_distribution<double> nd(0.0, 1.0);
        Eigen::MatrixXcd A(n, n); for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) A(i, j) = comp(nd(g), nd(g));
        Eigen::HouseholderQR<Eigen::MatrixXcd> qr(A); Eigen::MatrixXcd U = qr.householderQ() * Eigen::MatrixXcd::Identity(n, n);
        Eigen::MatrixXcd expect = U.adjoint() * o1.Q * U;
        std::cout << pr.first << " " << pr.second << " Nproj " << n << "/" << o2.Nproj << "  |Q'-U^dag Q U|/|Q| = " << (o2.Q - expect).norm() / expect.norm()
                  << "  Q' Hermitian dev " << (o2.Q - o2.Q.adjoint()).norm() / o2.Q.norm() << "  dln " << std::abs(o1.lnabsdet - o2.lnabsdet)
                  << "  sign " << o1.sign << "/" << o2.sign << "  sign_imag " << o2.sign_imag << "  max|Im Q'| " << o2.Q.imag().cwiseAbs().maxCoeff() << "\n";
    }
}
