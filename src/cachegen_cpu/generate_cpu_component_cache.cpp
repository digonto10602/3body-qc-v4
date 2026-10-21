// ============================================================================
// generate_cpu_component_cache.cpp
//
// CPU component cache generator for KKpi F3^{-1} ingredients.
// Built from the physics path in k3iso01_swave_cpu_cache_qc_omp_package,
// especially generate_swave_cpu_cache_qc.cpp and the v32f full-F3inv cache
// routines.  This version saves each ingredient into its own binary cache:
//
//   plm_config, klm_config, K2inv, F2, G, F3, F3inv, Vsel, projected_F3inv
//
// Default production physics follows the active KKpi defaults:
//   waves_vec_1={0,1}, waves_vec_2={0}, scatter1_00=4.04,
//   scatter1_10=-43.2, scatter2_00=4.12, xi=3.444, parity=-1.
// Set waves_vec_1=0 and waves_vec_2=0 in the config to reproduce the original
// s-wave-only CPU package behavior.
// ============================================================================

#define V32F_DISABLE_MINUIT 1
#define V31L_DISABLE_MINUIT 1
#define V31R_DISABLE_MINUIT 1
#define V31S_DISABLE_MINUIT 1
#define V31T_DISABLE_MINUIT 1

#include "K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;
using namespace k3df_fit_v32f;

static std::string trim_copy(std::string s) {
    while(!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
    while(!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
    return s;
}

static std::map<std::string,std::string> read_kv(const std::string& path) {
    std::ifstream in(path);
    if(!in) throw std::runtime_error("Could not open config: " + path);
    std::map<std::string,std::string> kv;
    std::string line;
    while(std::getline(in,line)) {
        const auto hash = line.find('#');
        if(hash != std::string::npos) line = line.substr(0, hash);
        const auto eq = line.find('=');
        if(eq == std::string::npos) continue;
        std::string k = trim_copy(line.substr(0,eq));
        std::string v = trim_copy(line.substr(eq+1));
        if(!k.empty()) kv[k] = v;
    }
    return kv;
}

static std::string gs(const std::map<std::string,std::string>& kv, const std::string& k, const std::string& d) {
    auto it=kv.find(k); return it==kv.end()?d:it->second;
}
static double gd(const std::map<std::string,std::string>& kv, const std::string& k, double d) {
    auto it=kv.find(k); return it==kv.end()?d:std::stod(it->second);
}
static int gi(const std::map<std::string,std::string>& kv, const std::string& k, int d) {
    auto it=kv.find(k); return it==kv.end()?d:std::stoi(it->second);
}

static std::vector<std::string> split_ws(std::string s) {
    for(char& c:s) if(c==',' || c=='{' || c=='}') c=' ';
    std::istringstream is(s);
    std::vector<std::string> out; std::string x;
    while(is>>x) out.push_back(x);
    return out;
}
static std::vector<double> split_doubles(const std::string& s, const std::vector<double>& dflt) {
    std::vector<std::string> toks = split_ws(s);
    if(toks.empty()) return dflt;
    std::vector<double> out; out.reserve(toks.size());
    for(const auto& t: toks) out.push_back(std::stod(t));
    return out;
}
static std::vector<int> split_ints(const std::string& s, const std::vector<int>& dflt) {
    std::vector<std::string> toks = split_ws(s);
    if(toks.empty()) return dflt;
    std::vector<int> out; out.reserve(toks.size());
    for(const auto& t: toks) out.push_back(std::stoi(t));
    return out;
}
static std::string ltag(double L) {
    double r = std::round(L);
    if(std::abs(L-r)<1e-9) return std::to_string(static_cast<int>(r));
    std::ostringstream os; os << std::fixed << std::setprecision(6) << L;
    std::string s=os.str();
    while(s.size()>1 && s.back()=='0') s.pop_back();
    if(!s.empty() && s.back()=='.') s.pop_back();
    for(char& c:s) if(c=='.') c='p';
    return s;
}
static std::string now_utc() {
    auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream os; os << std::put_time(std::gmtime(&t), "%Y-%m-%dT%H:%M:%SZ");
    return os.str();
}

struct Config {
    std::vector<double> Lbyas_values = {20.0, 24.0};
    std::vector<std::string> irreps = {"000_A1m", "100_A2", "110_A2", "111_A2", "200_A2"};
    double xi = 3.444;
    int coarseN = 1000;
    double Ecm_min = 0.26310;
    double Ecm_max = 0.36000;
    int omp_threads = 8;
    int print_every_percent = 10;
    std::string output_root = "output/cpu_component_cache";
    std::vector<int> waves_vec_1 = {0,1};
    std::vector<int> waves_vec_2 = {0};
    FitSettings settings;
};

static Config config_from_kv(const std::map<std::string,std::string>& kv) {
    Config cfg;
    cfg.Lbyas_values = split_doubles(gs(kv,"Lbyas_values","20 24"), cfg.Lbyas_values);
    cfg.irreps = split_ws(gs(kv,"irreps","000_A1m 100_A2 110_A2 111_A2 200_A2"));
    cfg.xi = gd(kv,"xi", gd(kv,"xival", cfg.xi));
    cfg.coarseN = gi(kv,"coarseN", cfg.coarseN);
    cfg.Ecm_min = gd(kv,"Ecm_min", cfg.Ecm_min);
    cfg.Ecm_max = gd(kv,"Ecm_max", cfg.Ecm_max);
    cfg.omp_threads = gi(kv,"omp_threads", cfg.omp_threads);
    cfg.print_every_percent = gi(kv,"print_every_percent", cfg.print_every_percent);
    cfg.output_root = gs(kv,"output_root", cfg.output_root);
    cfg.waves_vec_1 = split_ints(gs(kv,"waves_vec_1","0 1"), cfg.waves_vec_1);
    cfg.waves_vec_2 = split_ints(gs(kv,"waves_vec_2","0"), cfg.waves_vec_2);

    FitSettings s;
    s.xival = cfg.xi;
    s.coarseN = cfg.coarseN;
    s.omp_threads = cfg.omp_threads;
    s.debug = gs(kv,"debug","n").empty() ? 'n' : gs(kv,"debug","n")[0];
    s.atmpi = gd(kv,"atmpi", s.atmpi);
    s.atmK = gd(kv,"atmK", s.atmK);
    s.eta_1 = gd(kv,"eta_1", s.eta_1);
    s.eta_2 = gd(kv,"eta_2", s.eta_2);
    s.alpha = gd(kv,"alpha", s.alpha);
    s.epsilon_h = gd(kv,"epsilon_h", s.epsilon_h);
    s.max_shell_num = gd(kv,"max_shell_num", s.max_shell_num);
    s.tolerance = gd(kv,"tolerance", s.tolerance);
    s.parity = gi(kv,"parity", s.parity);
    s.eig_tol = gd(kv,"eig_tol", s.eig_tol);
    s.norm_tol = gd(kv,"norm_tol", s.norm_tol);
    s.proj_tol = gd(kv,"proj_tol", s.proj_tol);
    s.Q0norm = gi(kv,"Q0norm", s.Q0norm ? 1 : 0) != 0;
    s.sort_orbit_flag = gi(kv,"sort_orbit_flag", s.sort_orbit_flag ? 1 : 0) != 0;
    for(int a=0; a<4; ++a) for(int b=0; b<3; ++b) {
        s.scatter_params_1[a][b] = gd(kv,"scatter1_"+std::to_string(a)+std::to_string(b),s.scatter_params_1[a][b]);
        s.scatter_params_2[a][b] = gd(kv,"scatter2_"+std::to_string(a)+std::to_string(b),s.scatter_params_2[a][b]);
    }
    s.waves_vec_1 = cfg.waves_vec_1;
    s.waves_vec_2 = cfg.waves_vec_2;
    cfg.settings = s;

    if(std::find(cfg.waves_vec_1.begin(), cfg.waves_vec_1.end(), 1) != cfg.waves_vec_1.end()) {
        if(!std::isfinite(s.scatter_params_1[1][0]) || std::abs(s.scatter_params_1[1][0]) < 1.0e-300) {
            throw std::runtime_error("p-wave ell=1 is active but scatter1_10 is zero/non-finite");
        }
    }
    return cfg;
}

static void write_exact(std::ofstream& os, const void* p, std::size_t n) {
    os.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n));
    if(!os) throw std::runtime_error("binary write failed");
}
template<class T> static void write_scalar(std::ofstream& os, const T& x) { write_exact(os,&x,sizeof(T)); }
static void write_string(std::ofstream& os, const std::string& s) {
    std::uint32_t n = static_cast<std::uint32_t>(s.size());
    write_scalar(os,n);
    if(n) write_exact(os,s.data(),n);
}
static void write_matrix_record(std::ofstream& os, int row, double Ecm, const Eigen::MatrixXcd& M) {
    write_scalar(os, static_cast<std::int32_t>(row));
    write_scalar(os, Ecm);
    write_scalar(os, static_cast<std::int32_t>(M.rows()));
    write_scalar(os, static_cast<std::int32_t>(M.cols()));
    for(int r=0; r<M.rows(); ++r) for(int c=0; c<M.cols(); ++c) {
        const double re=M(r,c).real(); const double im=M(r,c).imag();
        write_scalar(os,re); write_scalar(os,im);
    }
}
static Eigen::MatrixXcd config_to_matrix(const std::vector<std::vector<comp>>& cfg) {
    const int cols = cfg.empty() ? 0 : static_cast<int>(cfg[0].size());
    int rows = static_cast<int>(cfg.size());
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(rows, cols);
    for(int r=0; r<rows; ++r) for(int c=0; c<static_cast<int>(cfg[std::size_t(r)].size()); ++c) M(r,c)=cfg[std::size_t(r)][std::size_t(c)];
    return M;
}

struct ComponentEntry {
    int row=-1;
    double Ecm=NAN;
    double En=NAN;
    int success=0;
    std::string error="OK";
    int dim1=0, dim2=0, total_dim=0, proj_dim=0;
    Eigen::MatrixXcd plm_config, klm_config, K2inv, F2, G, F3, F3inv, Vsel, projected_F3inv;
};

static ComponentEntry build_component_entry(int row, double Ecm, const MomentumIrrepSpec& spec, const FitSettings& s, const PhysicsParams& par) {
    ComponentEntry out; out.row=row; out.Ecm=Ecm;
    try {
        const std::vector<int> nnP_vec = {spec.nnP[0], spec.nnP[1], spec.nnP[2]};
        const comp pi = std::acos(-1.0);
        const double L = par.L();
        const comp twopibyL = comp(2.0,0.0)*pi/comp(L,0.0);
        std::vector<comp> total_P(3), nnP_config(3);
        for(int a=0; a<3; ++a) { total_P[a]=twopibyL*double(nnP_vec[a]); nnP_config[a]=comp(nnP_vec[a],0.0); }
        const comp Ecm_c(Ecm,0.0);
        const comp En_c = Ecm_to_E(Ecm_c,total_P);
        out.En = En_c.real();

        std::vector<std::vector<comp>> plm_config(5), klm_config(5);
        std::vector<std::vector<int>> np_config(5), nk_config(5);
        config_maker_4_momentum_first(plm_config,np_config,par.waves_vec_1,En_c,total_P,par.atmK,par.atmK,par.atmpi,L,par.epsilon_h,par.max_shell_num,par.tolerance);
        config_maker_4_momentum_first(klm_config,nk_config,par.waves_vec_2,En_c,total_P,par.atmpi,par.atmK,par.atmK,L,par.epsilon_h,par.max_shell_num,par.tolerance);
        out.plm_config = config_to_matrix(plm_config);
        out.klm_config = config_to_matrix(klm_config);
        out.dim1 = int(plm_config[0].size());
        out.dim2 = int(klm_config[0].size());
        out.total_dim = out.dim1 + out.dim2;
        if(out.total_dim <= 0) { out.error="EMPTY_CONFIG"; return out; }

        out.F2.resize(out.total_dim,out.total_dim);
        out.G.resize(out.total_dim,out.total_dim);
        out.K2inv.resize(out.total_dim,out.total_dim);
        F2_2plus1_mat(out.F2,En_c,plm_config,klm_config,total_P,par.atmK,par.atmpi,L,par.alpha,par.epsilon_h,par.max_shell_num,par.Q0norm);
        K2inv_EREord2_2plus1_mat(out.K2inv,par.eta_1,par.eta_2,par.scatter_params_1,par.scatter_params_2,En_c,plm_config,klm_config,total_P,par.atmK,par.atmpi,par.epsilon_h,L);
        G_2plus1_mat(out.G,En_c,plm_config,klm_config,total_P,par.atmK,par.atmpi,L,par.alpha,par.epsilon_h,par.max_shell_num,par.Q0norm);

        const Eigen::MatrixXcd H = out.K2inv + out.F2 + out.G;
        Eigen::PartialPivLU<Eigen::MatrixXcd> luH(H);
        Eigen::MatrixXcd X1 = luH.solve(out.F2);
        out.F3 = (out.F2/comp(3.0,0.0)) - out.F2*X1;

        CachedProjectorV30q pcache = get_projector_cached_v30q(plm_config,np_config,klm_config,nk_config,spec.irrep,nnP_config,par);
        out.Vsel = pcache.Vsel;
        out.proj_dim = int(out.Vsel.cols());
        if(out.proj_dim <= 0) { out.error="ZERO_PROJECTED_DIM"; return out; }

        double f3_scale = 0.0;
        for(int r=0; r<out.F3.rows(); ++r) for(int c=0; c<out.F3.cols(); ++c) {
            const double a = std::abs(out.F3(r,c));
            if(std::isfinite(a) && a > f3_scale) f3_scale = a;
        }
        if(!(f3_scale > 0.0) || !std::isfinite(f3_scale)) { out.error="BAD_F3_SCALE"; return out; }
        Eigen::MatrixXcd F3_scaled = out.F3 / f3_scale;
        Eigen::PartialPivLU<Eigen::MatrixXcd> luF3(F3_scaled);
        out.F3inv = luF3.solve(Eigen::MatrixXcd::Identity(out.total_dim,out.total_dim));
        out.F3inv /= f3_scale;
        out.projected_F3inv = out.Vsel.adjoint() * out.F3inv * out.Vsel;
        if(!out.K2inv.allFinite() || !out.F2.allFinite() || !out.G.allFinite() || !out.F3.allFinite() || !out.F3inv.allFinite() || !out.Vsel.allFinite() || !out.projected_F3inv.allFinite()) {
            out.error="NONFINITE_COMPONENT"; return out;
        }
        out.success=1; out.error="OK";
    } catch(const std::exception& e) {
        out.success=0; out.error=e.what();
    }
    return out;
}

static std::string object_path(const std::string& outdir, double L, double xi, const std::string& irrep, int coarseN, const std::string& object) {
    std::ostringstream os;
    os << outdir << "/cpu_Lbyas" << ltag(L) << "_xi" << ltag(xi) << "_irrep" << irrep << "_coarse" << coarseN << "_" << object << ".bin";
    return os.str();
}

struct Writers {
    std::map<std::string,std::string> path;
    std::map<std::string,std::ofstream> os;
};
static const std::vector<std::string>& object_names() {
    static const std::vector<std::string> names = {"plm_config","klm_config","K2inv","F2","G","F3","F3inv","Vsel","projected_F3inv"};
    return names;
}
static Writers open_writers(const std::string& outdir, double L, double xi, const std::string& irrep, int coarseN, std::size_t row_count) {
    fs::create_directories(outdir);
    Writers w;
    for(const auto& name: object_names()) {
        std::string p = object_path(outdir,L,xi,irrep,coarseN,name);
        w.path[name]=p;
        auto& os = w.os[name];
        os.open(p, std::ios::binary | std::ios::trunc);
        if(!os) throw std::runtime_error("Could not open component cache: "+p);
        os.write("KKPI_COMPONENT_CACHE_V001",25);
        write_string(os,name);
        write_string(os,"cpu_v33d_component_cachegen");
        write_scalar(os,L);
        write_scalar(os,xi);
        write_string(os,irrep);
        write_scalar(os,static_cast<std::int32_t>(coarseN));
        write_scalar(os,static_cast<std::int32_t>(row_count));
    }
    return w;
}
static void write_meta(const std::string& path, const std::string& object, const Config& cfg, const std::string& irrep, double L, long long attempted, long long ok, double sec, const std::string& started, const std::string& finished) {
    std::ofstream m(path + ".meta.json");
    if(!m) throw std::runtime_error("Could not write meta: "+path+".meta.json");
    m << std::setprecision(17);
    m << "{\n";
    m << "  \"version\": \"v33d_cpu_component_cachegen_all_ingredients\",\n";
    m << "  \"backend\": \"cpu\",\n";
    m << "  \"component_cache_file_format\": \"KKPI_COMPONENT_CACHE_V001\",\n";
    m << "  \"cache_kind\": \"" << object << "\",\n";
    m << "  \"cache_file\": \"" << path << "\",\n";
    m << "  \"irrep_label\": \"" << irrep << "\",\n";
    m << "  \"Lbyas\": " << L << ",\n";
    m << "  \"xi\": " << cfg.xi << ",\n";
    m << "  \"spatial_L\": " << L*cfg.xi << ",\n";
    m << "  \"Ecm_min\": " << cfg.Ecm_min << ",\n";
    m << "  \"Ecm_max\": " << cfg.Ecm_max << ",\n";
    m << "  \"coarseN\": " << cfg.coarseN << ",\n";
    m << "  \"waves_vec_1\": ["; for(std::size_t i=0;i<cfg.waves_vec_1.size();++i){ if(i)m<<", "; m<<cfg.waves_vec_1[i]; } m << "],\n";
    m << "  \"waves_vec_2\": ["; for(std::size_t i=0;i<cfg.waves_vec_2.size();++i){ if(i)m<<", "; m<<cfg.waves_vec_2[i]; } m << "],\n";
    m << "  \"scatter_params\": {\"scatter1_00\": " << cfg.settings.scatter_params_1[0][0] << ", \"scatter1_10\": " << cfg.settings.scatter_params_1[1][0] << ", \"scatter2_00\": " << cfg.settings.scatter_params_2[0][0] << "},\n";
    m << "  \"attempted_rows\": " << attempted << ",\n";
    m << "  \"success_rows\": " << ok << ",\n";
    m << "  \"elapsed_sec\": " << sec << ",\n";
    m << "  \"started_utc\": \"" << started << "\",\n";
    m << "  \"finished_utc\": \"" << finished << "\"\n";
    m << "}\n";
}

int main(int argc, char** argv) {
    if(argc != 2) {
        std::cerr << "Usage:\n  " << argv[0] << " configs/config_cpu_component_cachegen.in\n";
        return 2;
    }
    try {
        const std::string started = now_utc();
        const auto t0 = std::chrono::steady_clock::now();
        Config cfg = config_from_kv(read_kv(argv[1]));
        if(cfg.coarseN <= 1) throw std::runtime_error("coarseN must be > 1");
        if(!(cfg.Ecm_max > cfg.Ecm_min)) throw std::runtime_error("Ecm_max must be > Ecm_min");
#ifdef _OPENMP
        if(cfg.omp_threads > 0) omp_set_num_threads(cfg.omp_threads);
#endif
        std::cout << std::setprecision(17);
        std::cout << "[cpu-component-cachegen] config=" << argv[1] << "\n";
        std::cout << "[cpu-component-cachegen] output_root=" << cfg.output_root << "\n";
        std::cout << "[cpu-component-cachegen] coarseN=" << cfg.coarseN << " Ecm=[" << cfg.Ecm_min << "," << cfg.Ecm_max << "]\n";
        std::cout << "[cpu-component-cachegen] waves_vec_1={"; for(std::size_t i=0;i<cfg.waves_vec_1.size();++i){ if(i) std::cout<<","; std::cout<<cfg.waves_vec_1[i]; } std::cout << "} waves_vec_2={"; for(std::size_t i=0;i<cfg.waves_vec_2.size();++i){ if(i) std::cout<<","; std::cout<<cfg.waves_vec_2[i]; } std::cout << "}\n";
        fs::create_directories(cfg.output_root);
        const double del = (cfg.Ecm_max - cfg.Ecm_min) / double(cfg.coarseN - 1);
        long long grand_attempted=0, grand_ok=0;

        for(double L: cfg.Lbyas_values) {
            FitSettings sL = cfg.settings;
            sL.Lval = L;
            sL.xival = cfg.xi;
            sL.coarseN = cfg.coarseN;
            sL.waves_vec_1 = cfg.waves_vec_1;
            sL.waves_vec_2 = cfg.waves_vec_2;
            PhysicsParams par = make_base_physics(sL);
            par.waves_vec_1 = cfg.waves_vec_1;
            par.waves_vec_2 = cfg.waves_vec_2;
            const std::string outdir = cfg.output_root + "/cpu_Lbyas" + ltag(L) + "_component_cache";
            fs::create_directories(outdir);
            for(const auto& irrep: cfg.irreps) {
                const auto t_ir0 = std::chrono::steady_clock::now();
                MomentumIrrepSpec spec = parse_label(irrep);
                Writers writers = open_writers(outdir,L,cfg.xi,irrep,cfg.coarseN,static_cast<std::size_t>(cfg.coarseN));
                long long ok=0;
                std::cout << "[irrep-start] L=" << L << " irrep=" << irrep << " outdir=" << outdir << "\n";
                int next_pct = cfg.print_every_percent;
                for(int i=0;i<cfg.coarseN;++i) {
                    const double Ecm = cfg.Ecm_min + double(i)*del;
                    ComponentEntry e = build_component_entry(i,Ecm,spec,sL,par);
                    ++grand_attempted;
                    if(e.success) {
                        ++ok; ++grand_ok;
                        write_matrix_record(writers.os["plm_config"], i, Ecm, e.plm_config);
                        write_matrix_record(writers.os["klm_config"], i, Ecm, e.klm_config);
                        write_matrix_record(writers.os["K2inv"], i, Ecm, e.K2inv);
                        write_matrix_record(writers.os["F2"], i, Ecm, e.F2);
                        write_matrix_record(writers.os["G"], i, Ecm, e.G);
                        write_matrix_record(writers.os["F3"], i, Ecm, e.F3);
                        write_matrix_record(writers.os["F3inv"], i, Ecm, e.F3inv);
                        write_matrix_record(writers.os["Vsel"], i, Ecm, e.Vsel);
                        write_matrix_record(writers.os["projected_F3inv"], i, Ecm, e.projected_F3inv);
                    } else {
                        std::cout << "[row-fail] L=" << L << " irrep=" << irrep << " i=" << i << " Ecm=" << Ecm << " error=" << e.error << "\n";
                    }
                    int pct = int(100.0 * double(i+1) / double(cfg.coarseN));
                    if(pct >= next_pct) {
                        std::cout << "[progress] L=" << L << " irrep=" << irrep << " " << pct << "% ok=" << ok << "/" << (i+1) << "\n";
                        next_pct += cfg.print_every_percent;
                    }
                }
                for(auto& kv: writers.os) kv.second.close();
                const auto t_ir1 = std::chrono::steady_clock::now();
                const double sec = std::chrono::duration<double>(t_ir1-t_ir0).count();
                const std::string finished_ir = now_utc();
                for(const auto& name: object_names()) write_meta(writers.path[name],name,cfg,irrep,L,cfg.coarseN,ok,sec,started,finished_ir);
                std::cout << "[irrep-done] L=" << L << " irrep=" << irrep << " ok_rows=" << ok << "/" << cfg.coarseN << " sec=" << sec << "\n";
            }
        }
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        std::cout << "[DONE] CPU component cachegen attempted_rows=" << grand_attempted << " success_rows=" << grand_ok << " elapsed_sec=" << sec << " output_root=" << cfg.output_root << "\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "[cpu-component-cachegen ERROR] " << e.what() << "\n";
        return 1;
    }
}
