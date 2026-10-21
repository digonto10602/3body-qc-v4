// generate_swave_cpu_cache_qc.cpp
// CPU-only s-wave F3inv/Vsel cache + QC diagnostic generator.
//
// Requested behavior:
//   - waves_vec_1 = {0}, waves_vec_2 = {0}
//   - Kdf3 uses K3iso0 and K3iso1 only; K3B=0, K3E=0
//   - coarseN default 2000
//   - all default irreps: 000_A1m, 100_A2, 110_A2, 111_A2, 200_A2
//   - Lbyas volumes: 20, 24
//   - matrix generation on CPU
//   - save F3inv/Vsel caches for all L/irreps in one binary file: swave_L_20_24_cache.bin
//   - write text data with QC_isotropic and det(scaled_projected_QC)

#define V32F_DISABLE_MINUIT 1
#define V31L_DISABLE_MINUIT 1
#define V31R_DISABLE_MINUIT 1
#define V31S_DISABLE_MINUIT 1
#define V31T_DISABLE_MINUIT 1

#include "K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <atomic>
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
#include <limits>
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
static std::string swave_sanitize_error(std::string e) {
    if(e.empty()) return "OK";
    for(char& c:e) if(std::isspace(static_cast<unsigned char>(c))) c='_';
    return e;
}
static std::string ltag(double L) {
    double r = std::round(L);
    if(std::abs(L-r)<1e-9) return std::to_string(static_cast<int>(r));
    std::ostringstream os; os << std::fixed << std::setprecision(3) << L;
    std::string s=os.str(); for(char& c:s) if(c=='.') c='p'; return s;
}

struct Config {
    std::vector<double> Lbyas_values = {20.0, 24.0};
    std::vector<std::string> irreps = {"000_A1m", "100_A2", "110_A2", "111_A2", "200_A2"};
    double xi = 3.444;
    int coarseN = 2000;
    double Ecm_min = 0.26301;
    double Ecm_max = 0.36;
    int omp_threads = 8;
    double K3iso0 = 73735.840894011912;
    double K3iso1 = -972421.14060757787;
    double K3B = 0.0;
    double K3E = 0.0;
    double const_norm_power = 6.0;
    std::string cache_file = "output/cache/swave_L_20_24_cache.bin";
    std::string output_dir = "output/qc_swave";
    std::string combined_output_file = "output/qc_swave/swave_L_20_24_QC_isotropic_det_scaled_projected_QC.dat";
    std::string meta_file = "";
    int force_rebuild = 1;
    int print_every_percent = 10;
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
    cfg.K3iso0 = gd(kv,"K3iso0", cfg.K3iso0);
    cfg.K3iso1 = gd(kv,"K3iso1", cfg.K3iso1);
    cfg.K3B = 0.0;
    cfg.K3E = 0.0;
    cfg.const_norm_power = gd(kv,"const_norm_power", cfg.const_norm_power);
    cfg.cache_file = gs(kv,"cache_file", cfg.cache_file);
    cfg.output_dir = gs(kv,"output_dir", cfg.output_dir);
    cfg.combined_output_file = gs(kv,"combined_output_file", cfg.combined_output_file);
    cfg.meta_file = gs(kv,"meta_file", cfg.cache_file + ".meta.json");
    cfg.force_rebuild = gi(kv,"force_rebuild", cfg.force_rebuild);
    cfg.print_every_percent = gi(kv,"print_every_percent", cfg.print_every_percent);

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
    // Hard request: s-wave only.
    s.waves_vec_1 = {0};
    s.waves_vec_2 = {0};
    cfg.settings = s;
    return cfg;
}

static void write_bytes(std::ofstream& os, const void* p, std::size_t n) {
    os.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n));
    if(!os) throw std::runtime_error("binary cache write failed");
}
template<class T> static void write_scalar(std::ofstream& os, const T& x) { write_bytes(os,&x,sizeof(T)); }
static void write_string_bin(std::ofstream& os, const std::string& s) {
    std::uint32_t n = static_cast<std::uint32_t>(s.size());
    write_scalar(os,n);
    if(n) write_bytes(os,s.data(),n);
}
static void write_matrix_bin(std::ofstream& os, const Eigen::MatrixXcd& M) {
    std::int32_t rows = static_cast<std::int32_t>(M.rows());
    std::int32_t cols = static_cast<std::int32_t>(M.cols());
    write_scalar(os,rows); write_scalar(os,cols);
    for(int c=0; c<M.cols(); ++c) for(int r=0; r<M.rows(); ++r) {
        double re = M(r,c).real(); double im = M(r,c).imag();
        write_scalar(os,re); write_scalar(os,im);
    }
}

static double determinant_real(const Eigen::MatrixXcd& M) {
    if(M.rows()==0 || M.rows()!=M.cols()) return std::numeric_limits<double>::quiet_NaN();
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu(M);
    return lu.determinant().real();
}

static double signed_logabs_det(const Eigen::MatrixXcd& M, double* logabs_out=nullptr, double* sign_out=nullptr) {
    if(M.rows()==0 || M.rows()!=M.cols()) {
        if(logabs_out) *logabs_out = NAN;
        if(sign_out) *sign_out = NAN;
        return NAN;
    }
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu(M);
    comp detc = lu.determinant();
    double sign = (detc.real()>0.0) ? 1.0 : ((detc.real()<0.0) ? -1.0 : 0.0);
    const auto& LU = lu.matrixLU();
    double lga = 0.0;
    bool ok = LU.rows()==LU.cols() && LU.rows()>0;
    for(int d=0; d<LU.rows(); ++d) {
        double a = std::abs(LU(d,d));
        if(!(a>0.0) || !std::isfinite(a)) { ok=false; break; }
        lga += std::log(a);
    }
    if(!ok) {
        double ad = std::abs(detc);
        lga = (ad>0.0 && std::isfinite(ad)) ? std::log(ad) : std::numeric_limits<double>::quiet_NaN();
    }
    if(logabs_out) *logabs_out = lga;
    if(sign_out) *sign_out = sign;
    return (std::isfinite(lga) && sign != 0.0) ? sign*lga : std::numeric_limits<double>::quiet_NaN();
}

static comp k3iso_scalar_from_Ecm(double Ecm, const Config& cfg, const PhysicsParams& par) {
    const double Mtot = 2.0*par.atmK + par.atmpi;
    const comp Delta = comp((Ecm*Ecm)/(Mtot*Mtot) - 1.0, 0.0);
    return comp(cfg.K3iso0,0.0) + comp(cfg.K3iso1,0.0)*Delta;
}

static comp compute_F3iso_inv_from_F3inv_full(const ProjectedQCCacheEntry& e) {
    const int A = int(e.plm_config[0].size());
    const int B = int(e.klm_config[0].size());
    const int N = A + B;
    if(N <= 0 || e.F3inv_full.rows() != N || e.F3inv_full.cols() != N) return comp(NAN,NAN);
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu(e.F3inv_full);
    Eigen::MatrixXcd F3 = lu.solve(Eigen::MatrixXcd::Identity(N,N));
    Eigen::VectorXcd v_iso(N);
    for(int r=0; r<A; ++r) v_iso(r) = comp(1.0,0.0);
    for(int r=0; r<B; ++r) v_iso(A+r) = comp(1.0/std::sqrt(2.0),0.0);
    comp F3iso = (v_iso.transpose() * F3 * v_iso)(0,0);
    if(std::abs(F3iso) == 0.0) return comp(NAN,NAN);
    return comp(1.0,0.0) / F3iso;
}

struct RowDiagnostics {
    int A=0, B=0;
    double Delta=NAN;
    double F3isoInv_re=NAN;
    double K3iso_re=NAN;
    double QCiso_re=NAN;
    double det_scaled=NAN;
    double signed_lg=NAN;
    double logabs=NAN;
    double det_sign=NAN;
    double closest_eig=NAN;
    double min_abs_eig=NAN;
    double min_sv=NAN;
};

struct RowResult {
    int row = -1;
    double Ecm = std::numeric_limits<double>::quiet_NaN();
    ProjectedQCCacheEntry e;
    RowDiagnostics d;
};

static RowDiagnostics compute_diagnostics(const ProjectedQCCacheEntry& e, const Config& cfg, const PhysicsParams& par, double cscale) {
    RowDiagnostics d;
    d.A = int(e.plm_config.empty()?0:e.plm_config[0].size());
    d.B = int(e.klm_config.empty()?0:e.klm_config[0].size());
    const double Mtot = 2.0*par.atmK + par.atmpi;
    d.Delta = (e.Ecm*e.Ecm)/(Mtot*Mtot) - 1.0;
    if(!e.success) return d;

    comp F3isoInv = compute_F3iso_inv_from_F3inv_full(e);
    comp Kiso = k3iso_scalar_from_Ecm(e.Ecm,cfg,par);
    comp QCiso = F3isoInv + Kiso;
    d.F3isoInv_re = F3isoInv.real();
    d.K3iso_re = Kiso.real();
    d.QCiso_re = QCiso.real();

    Eigen::MatrixXcd K0_proj = make_K3_projected(e.plm_config,e.klm_config,e.total_P,e.En_c,par,e.Vsel,1.0,0.0,0.0,0.0,'n');
    Eigen::MatrixXcd K1_proj = make_K3_projected(e.plm_config,e.klm_config,e.total_P,e.En_c,par,e.Vsel,0.0,1.0,0.0,0.0,'n');
    Eigen::MatrixXcd Q = e.F3inv_proj + comp(cfg.K3iso0,0.0)*K0_proj + comp(cfg.K3iso1,0.0)*K1_proj;
    Eigen::MatrixXcd Q_scaled = Q / comp(cscale,0.0);
    d.det_scaled = determinant_real(Q_scaled);
    d.signed_lg = signed_logabs_det(Q_scaled, &d.logabs, &d.det_sign);

    double herm=NAN; int eigidx=-1; double eig=NAN;
    if(closest_zero_eigenvalue(Q_scaled,eig,eigidx,herm)) {
        d.closest_eig = eig;
        d.min_abs_eig = std::abs(eig);
    }
    try {
        Eigen::BDCSVD<Eigen::MatrixXcd> svd(Q_scaled, Eigen::ComputeThinU | Eigen::ComputeThinV);
        if(svd.singularValues().size() > 0) d.min_sv = svd.singularValues().minCoeff();
    } catch(...) {}
    return d;
}

static void write_meta_json(const Config& cfg, const std::string& started, const std::string& finished, long long total_rows, long long ok_rows, double elapsed_sec) {
    fs::create_directories(fs::path(cfg.meta_file).parent_path());
    std::ofstream m(cfg.meta_file);
    if(!m) throw std::runtime_error("Could not write meta: " + cfg.meta_file);
    m << std::setprecision(17);
    m << "{\n";
    m << "  \"version\": \"swave_cpu_cache_qc_v1\",\n";
    m << "  \"cache_kind\": \"F3inv_full_plus_Vsel_swave_cpu\",\n";
    m << "  \"cache_file\": \"" << cfg.cache_file << "\",\n";
    m << "  \"coarseN\": " << cfg.coarseN << ",\n";
    m << "  \"Ecm_min\": " << cfg.Ecm_min << ",\n";
    m << "  \"Ecm_max\": " << cfg.Ecm_max << ",\n";
    m << "  \"xi\": " << cfg.xi << ",\n";
    m << "  \"waves_vec_1\": [0],\n";
    m << "  \"waves_vec_2\": [0],\n";
    m << "  \"K3iso0\": " << cfg.K3iso0 << ",\n";
    m << "  \"K3iso1\": " << cfg.K3iso1 << ",\n";
    m << "  \"K3B\": 0,\n";
    m << "  \"K3E\": 0,\n";
    m << "  \"scatter1_00\": " << cfg.settings.scatter_params_1[0][0] << ",\n";
    m << "  \"scatter1_10\": " << cfg.settings.scatter_params_1[1][0] << ",\n";
    m << "  \"scatter2_00\": " << cfg.settings.scatter_params_2[0][0] << ",\n";
    m << "  \"Lbyas_values\": [";
    for(std::size_t i=0;i<cfg.Lbyas_values.size();++i){ if(i) m << ", "; m << cfg.Lbyas_values[i]; }
    m << "],\n";
    m << "  \"irreps\": [";
    for(std::size_t i=0;i<cfg.irreps.size();++i){ if(i) m << ", "; m << "\"" << cfg.irreps[i] << "\""; }
    m << "],\n";
    m << "  \"total_rows\": " << total_rows << ",\n";
    m << "  \"success_rows\": " << ok_rows << ",\n";
    m << "  \"started_utc_note\": \"" << started << "\",\n";
    m << "  \"finished_utc_note\": \"" << finished << "\",\n";
    m << "  \"elapsed_sec\": " << elapsed_sec << "\n";
    m << "}\n";
}

static std::string now_note() {
    auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream os; os << std::put_time(std::gmtime(&t), "%Y-%m-%dT%H:%M:%SZ");
    return os.str();
}

int main(int argc, char** argv) {
    if(argc != 2) {
        std::cerr << "Usage:\n  " << argv[0] << " configs/config_swave_cpu_cache_qc.in\n";
        return 1;
    }
    try {
        const std::string started = now_note();
        auto t0 = std::chrono::steady_clock::now();
        auto kv = read_kv(argv[1]);
        Config cfg = config_from_kv(kv);

        if(cfg.coarseN <= 1) throw std::runtime_error("coarseN must be > 1");
        if(!(cfg.Ecm_max > cfg.Ecm_min)) throw std::runtime_error("Ecm_max must be > Ecm_min");
        fs::create_directories(fs::path(cfg.cache_file).parent_path());
        fs::create_directories(cfg.output_dir);
        fs::create_directories(fs::path(cfg.combined_output_file).parent_path());

#ifdef _OPENMP
        if(cfg.omp_threads > 0) omp_set_num_threads(cfg.omp_threads);
#endif
        std::cout << std::setprecision(17);
        std::cout << "[swave-cpu-qc] config=" << argv[1] << "\n";
        std::cout << "[swave-cpu-qc] cache_file=" << cfg.cache_file << "\n";
        std::cout << "[swave-cpu-qc] combined_output_file=" << cfg.combined_output_file << "\n";
        std::cout << "[swave-cpu-qc] coarseN=" << cfg.coarseN << " Ecm=[" << cfg.Ecm_min << "," << cfg.Ecm_max << "]\n";
        std::cout << "[swave-cpu-qc] waves_vec_1={0} waves_vec_2={0}\n";
        std::cout << "[swave-cpu-qc] K3iso0=" << cfg.K3iso0 << " K3iso1=" << cfg.K3iso1 << " K3B=0 K3E=0\n";
        std::cout << "[swave-cpu-qc] L values:"; for(double L: cfg.Lbyas_values) std::cout << ' ' << L; std::cout << "\n";
        std::cout << "[swave-cpu-qc] irreps:"; for(const auto& x: cfg.irreps) std::cout << ' ' << x; std::cout << "\n";

        std::ofstream cache(cfg.cache_file, std::ios::binary | std::ios::trunc);
        if(!cache) throw std::runtime_error("Could not open binary cache for writing: " + cfg.cache_file);
        std::ofstream combined(cfg.combined_output_file);
        if(!combined) throw std::runtime_error("Could not open combined output: " + cfg.combined_output_file);
        combined << std::setprecision(17);
        combined << "# swave CPU-generated QC data\n";
        combined << "# cache_file " << cfg.cache_file << "\n";
        combined << "# waves_vec_1 0\n# waves_vec_2 0\n";
        combined << "# K3iso0 " << cfg.K3iso0 << "\n# K3iso1 " << cfg.K3iso1 << "\n# K3B 0\n# K3E 0\n";
        combined << "# columns: Lbyas irrep row Ecm success A B total_dim proj_dim Delta F3isoInv_re K3iso_re QC_isotropic_re det_scaled_projected_QC_re signed_logabsdet_scaled logabsdet_scaled det_sign closest_eig_scaled min_abs_eig_scaled min_sv_scaled error\n";

        const std::string header = "KKPI_SWAVE_L20_24_F3INV_VSEL_CPU_CACHE_V1\n";
        cache.write(header.data(), static_cast<std::streamsize>(header.size()));
        std::uint32_t version = 1;
        std::uint32_t nL = static_cast<std::uint32_t>(cfg.Lbyas_values.size());
        std::uint32_t nIr = static_cast<std::uint32_t>(cfg.irreps.size());
        write_scalar(cache, version);
        write_scalar(cache, nL);
        write_scalar(cache, nIr);
        write_scalar(cache, cfg.coarseN);
        write_scalar(cache, cfg.xi);
        write_scalar(cache, cfg.Ecm_min);
        write_scalar(cache, cfg.Ecm_max);
        write_scalar(cache, cfg.K3iso0);
        write_scalar(cache, cfg.K3iso1);
        double zero=0.0; write_scalar(cache, zero); write_scalar(cache, zero); // K3B, K3E
        for(double L: cfg.Lbyas_values) write_scalar(cache, L);
        for(const auto& ir: cfg.irreps) write_string_bin(cache, ir);

        long long total_rows = 0;
        long long ok_rows = 0;
        const long long total_expected = static_cast<long long>(cfg.Lbyas_values.size()) * static_cast<long long>(cfg.irreps.size()) * static_cast<long long>(cfg.coarseN);
        int global_next_pct = cfg.print_every_percent;

        for(std::size_t iL=0; iL<cfg.Lbyas_values.size(); ++iL) {
            const double Lbyas = cfg.Lbyas_values[iL];
            FitSettings sL = cfg.settings;
            sL.Lval = Lbyas;
            sL.xival = cfg.xi;
            sL.coarseN = cfg.coarseN;
            sL.waves_vec_1 = {0};
            sL.waves_vec_2 = {0};
            PhysicsParams par = make_base_physics(sL);
            par.waves_vec_1 = {0};
            par.waves_vec_2 = {0};
            const double cscale = std::pow(Lbyas * cfg.xi, cfg.const_norm_power);

            for(std::size_t ii=0; ii<cfg.irreps.size(); ++ii) {
                const std::string irrep = cfg.irreps[ii];
                MomentumIrrepSpec spec = parse_label(irrep);
                const std::string per_file = (fs::path(cfg.output_dir) / ("L" + ltag(Lbyas) + "_" + irrep + "_swave_QC.dat")).string();
                std::ofstream per(per_file);
                if(!per) throw std::runtime_error("Could not open per-irrep output: " + per_file);
                per << std::setprecision(17);
                per << "# swave CPU-generated QC data\n";
                per << "# Lbyas " << Lbyas << "\n# irrep " << irrep << "\n# coarseN " << cfg.coarseN << "\n";
                per << "# waves_vec_1 0\n# waves_vec_2 0\n";
                per << "# columns: row Ecm success A B total_dim proj_dim Delta F3isoInv_re K3iso_re QC_isotropic_re det_scaled_projected_QC_re signed_logabsdet_scaled logabsdet_scaled det_sign closest_eig_scaled min_abs_eig_scaled min_sv_scaled error\n";

                std::uint64_t block_magic = 0x5357415645424c4bULL; // SWAVEBLK
                write_scalar(cache, block_magic);
                write_scalar(cache, Lbyas);
                write_string_bin(cache, irrep);
                std::int32_t nrows_i32 = cfg.coarseN;
                write_scalar(cache, nrows_i32);
                write_scalar(cache, cscale);

                std::cout << "[block] L=" << Lbyas << " irrep=" << irrep << " cscale=" << cscale << " output=" << per_file << "\n";
#ifdef _OPENMP
                const int active_threads = (cfg.omp_threads > 0) ? cfg.omp_threads : omp_get_max_threads();
#else
                const int active_threads = 1;
#endif
                std::cout << "[block omp] L=" << Lbyas << " irrep=" << irrep
                          << " generating " << cfg.coarseN << " rows with " << active_threads
                          << " OpenMP thread(s); ordered file write happens after parallel generation.\n";

                std::vector<RowResult> rows(static_cast<std::size_t>(cfg.coarseN));
                std::atomic<int> block_done{0};
                int block_next_pct = cfg.print_every_percent;

#pragma omp parallel for schedule(dynamic) if(cfg.omp_threads != 1)
                for(int i=0; i<cfg.coarseN; ++i) {
                    RowResult rr;
                    rr.row = i;
                    const double t = double(i) / double(cfg.coarseN - 1);
                    const double Ecm = cfg.Ecm_min + t*(cfg.Ecm_max - cfg.Ecm_min);
                    rr.Ecm = Ecm;
                    try {
                        rr.e = build_cache_entry(i, Ecm, spec, sL, par, sL.debug);
                        rr.d = compute_diagnostics(rr.e, cfg, par, cscale);
                    } catch(const std::exception& ex) {
                        rr.e.label = irrep;
                        rr.e.spec = spec;
                        rr.e.i = i;
                        rr.e.Ecm = Ecm;
                        rr.e.success = 0;
                        rr.e.error = std::string("EXCEPTION_") + ex.what();
                        const double Mtot = 2.0*par.atmK + par.atmpi;
                        rr.d.Delta = (Ecm*Ecm)/(Mtot*Mtot) - 1.0;
                    } catch(...) {
                        rr.e.label = irrep;
                        rr.e.spec = spec;
                        rr.e.i = i;
                        rr.e.Ecm = Ecm;
                        rr.e.success = 0;
                        rr.e.error = "UNKNOWN_EXCEPTION";
                        const double Mtot = 2.0*par.atmK + par.atmpi;
                        rr.d.Delta = (Ecm*Ecm)/(Mtot*Mtot) - 1.0;
                    }
                    rows[static_cast<std::size_t>(i)] = std::move(rr);

                    const int done = ++block_done;
                    const int bpct = int(100.0 * double(done) / double(cfg.coarseN));
                    if(bpct >= block_next_pct) {
#pragma omp critical(swave_progress_print)
                        {
                            while(bpct >= block_next_pct) {
                                std::cout << "[block progress] L=" << Lbyas << " irrep=" << irrep
                                          << " " << block_next_pct << "% rows=" << done << "/" << cfg.coarseN << "\n";
                                block_next_pct += cfg.print_every_percent;
                            }
                        }
                    }
                }

                // File writes are deliberately single-threaded and ordered, so the one-file
                // binary cache is deterministic and free of write races.
                for(int i=0; i<cfg.coarseN; ++i) {
                    const RowResult& rr = rows[static_cast<std::size_t>(i)];
                    const ProjectedQCCacheEntry& e = rr.e;
                    const RowDiagnostics& d = rr.d;
                    const double Ecm = rr.Ecm;
                    if(e.success) ++ok_rows;
                    ++total_rows;

                    // Binary row record.
                    std::uint64_t row_magic = 0x5357415645524f57ULL; // SWAVEROW
                    write_scalar(cache, row_magic);
                    std::int32_t row_i32 = i;
                    std::int32_t success_i32 = e.success;
                    std::int32_t total_dim_i32 = e.total_dim;
                    std::int32_t proj_dim_i32 = e.proj_dim;
                    write_scalar(cache, row_i32);
                    write_scalar(cache, Ecm);
                    write_scalar(cache, e.En);
                    write_scalar(cache, success_i32);
                    write_scalar(cache, total_dim_i32);
                    write_scalar(cache, proj_dim_i32);
                    write_scalar(cache, d.QCiso_re);
                    write_scalar(cache, d.det_scaled);
                    write_scalar(cache, d.signed_lg);
                    write_scalar(cache, d.logabs);
                    write_scalar(cache, d.det_sign);
                    write_matrix_bin(cache, e.F3inv_full);
                    write_matrix_bin(cache, e.Vsel);
                    write_string_bin(cache, swave_sanitize_error(e.error));

                    per << i << ' ' << Ecm << ' ' << e.success << ' '
                        << d.A << ' ' << d.B << ' ' << e.total_dim << ' ' << e.proj_dim << ' '
                        << d.Delta << ' ' << d.F3isoInv_re << ' ' << d.K3iso_re << ' ' << d.QCiso_re << ' '
                        << d.det_scaled << ' ' << d.signed_lg << ' ' << d.logabs << ' ' << d.det_sign << ' '
                        << d.closest_eig << ' ' << d.min_abs_eig << ' ' << d.min_sv << ' '
                        << swave_sanitize_error(e.error) << "\n";
                    combined << Lbyas << ' ' << irrep << ' ' << i << ' ' << Ecm << ' ' << e.success << ' '
                        << d.A << ' ' << d.B << ' ' << e.total_dim << ' ' << e.proj_dim << ' '
                        << d.Delta << ' ' << d.F3isoInv_re << ' ' << d.K3iso_re << ' ' << d.QCiso_re << ' '
                        << d.det_scaled << ' ' << d.signed_lg << ' ' << d.logabs << ' ' << d.det_sign << ' '
                        << d.closest_eig << ' ' << d.min_abs_eig << ' ' << d.min_sv << ' '
                        << swave_sanitize_error(e.error) << "\n";

                    const int gpct = int(100.0 * double(total_rows) / double(total_expected));
                    if(gpct >= global_next_pct) {
                        std::cout << "[global progress] " << gpct << "% rows=" << total_rows << "/" << total_expected << "\n";
                        global_next_pct += cfg.print_every_percent;
                    }
                }
                per.close();
                std::cout << "[block done] L=" << Lbyas << " irrep=" << irrep
                          << " generated_parallel_rows=" << cfg.coarseN << " wrote_ordered_rows=" << cfg.coarseN << "\n";
            }
        }
        cache.close();
        combined.close();
        auto t1 = std::chrono::steady_clock::now();
        const double elapsed_sec = std::chrono::duration<double>(t1-t0).count();
        const std::string finished = now_note();
        write_meta_json(cfg, started, finished, total_rows, ok_rows, elapsed_sec);

        std::cout << "[swave-cpu-qc] wrote cache: " << cfg.cache_file << "\n";
        std::cout << "[swave-cpu-qc] wrote sidecar: " << cfg.meta_file << "\n";
        std::cout << "[swave-cpu-qc] wrote combined data: " << cfg.combined_output_file << "\n";
        std::cout << "[swave-cpu-qc] rows success=" << ok_rows << " total=" << total_rows << " elapsed_sec=" << elapsed_sec << "\n";
        std::cout << "[swave-cpu-qc] plot command:\n  python3 scripts/plot_swave_qc_allirreps.py " << cfg.combined_output_file << " --n-scale 100\n";
        std::cout << "[swave-cpu-qc] done\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "[swave-cpu-qc-error] " << e.what() << "\n";
        return 2;
    }
}
