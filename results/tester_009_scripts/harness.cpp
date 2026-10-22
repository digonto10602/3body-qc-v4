// tester harness (not production): exercises build_orbit_V self-checks and other irreps.  Mutants shadow fast_projector.hpp via -I order.
#include "fast_basis.hpp"
#include "fast_projector.hpp"
#include <iostream>
using namespace qc3_fast;
static double rel(const Eigen::MatrixXcd& a, const Eigen::MatrixXcd& b) { return (a - b).norm() / b.norm(); }
int main(int argc, char** argv) {
    std::string mode = argc > 1 ? argv[1] : "all";
    qc3_oracle::Input in;
    in.Ecm = 0.2915561; in.irrep = "100_A2";
    if (mode == "normal" || mode == "all") {
        for (auto pr : std::vector<std::pair<std::string, double>>{{"000_A1m", 0.34538}, {"100_A2", 0.2915561}, {"111_A2", 0.35103}, {"110_A2", 0.3400}, {"200_A2", 0.34738}}) {
            in.irrep = pr.first; in.Ecm = pr.second;
            auto k = fast_kinematics(in);
            try {
                Eigen::MatrixXd V = build_orbit_V(k);
                std::cout << "OK " << pr.first << " N=" << k.N << " Nproj=" << V.cols() << "\n";
            } catch (const std::exception& e) { std::cout << "THROW " << pr.first << ": " << e.what() << "\n"; }
        }
    }
    if (mode == "drop" || mode == "all") {
        in.irrep = "100_A2"; in.Ecm = 0.2915561;
        auto k = fast_kinematics(in);
        for (int row : {0, 5, k.N - 1, k.N / 2}) {
            auto k2 = k; k2.basis.erase(k2.basis.begin() + row);
            try { build_orbit_V(k2); std::cout << "NO THROW after dropping row " << row << "\n"; }
            catch (const std::exception& e) { std::cout << "THROW (drop row " << row << "): " << e.what() << "\n"; }
        }
    }
    if (mode == "irreps" || mode == "all") {
        // other irreps of the same groups (not production): orbit vs the existing eigen builder; irrep name overwritten after kinematics
        struct C { std::string base; double E; std::vector<std::string> irr; };
        for (auto c : std::vector<C>{{"100_A2", 0.31, {"A1", "A2", "B1", "B2", "E"}}, {"110_A2", 0.31, {"A1", "A2", "B1", "B2"}}, {"111_A2", 0.32, {"A1", "A2", "E"}},
                                      {"000_A1m", 0.32, {"A1g", "A2g", "A1u", "A2u", "Eg", "Eu", "T1g", "T1u", "T2g", "T2u"}}, {"200_A2", 0.31, {"A1", "A2", "B1", "B2", "E"}}}) {
            in.irrep = c.base; in.Ecm = c.E;
            auto k0 = fast_kinematics(in);
            for (auto ir : c.irr) {
                auto k = k0; k.spec.irrep = ir;
                try {
                    Eigen::MatrixXd Vo = build_orbit_V(k);
                    auto built = build_cached_projector_v30q(k.plm_config, k.np_config, k.klm_config, k.nk_config, k.spec.irrep, k.nnP_config, k.par);
                    Eigen::MatrixXcd Ve = built.Vsel, Vc = Vo.cast<comp>();
                    double d = (Ve.cols() == Vc.cols()) ? rel(Vc * Vc.adjoint(), Ve * Ve.adjoint()) : -1;
                    std::cout << "IRREP " << c.base << " as " << ir << " dim=" << fvproj_v30j::irrep_dimension(ir) << " N=" << k.N << " Nproj orbit=" << Vo.cols() << " eigen=" << Ve.cols() << " rel(VVdag)=" << d << "\n";
                } catch (const std::exception& e) { std::cout << "IRREP " << c.base << " as " << ir << " THROW: " << e.what() << "\n"; }
            }
        }
    }
    return 0;
}
