#include "K3df_minuit_fit_v32f_fullF3inv_QCfull_cached_classifier.hpp"
#include "v33g_runtime_k3basis_cache.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;
using namespace k3df_fit_v32f;
using comp = std::complex<double>;

namespace {

struct Args {
    fs::path gpu_cache_dir = "/media/digonto/Data/F3inv_cache";
    fs::path outdir = "output/oldscale_det_compare";
    fs::path plot_script = "scripts/plot_oldscale_det_compare.py";
    std::vector<double> Lbyas;
    std::vector<std::string> irreps;
    int points = 1000;
    int threads = 16;
    double xi = 3.444;
    bool k3df_zero = false;
    K3dfParameters k3{73735.840894095912, -972421.14060757787, 347174.05548116949, -1226756.7068845264};
};

struct GpuRowMeta {
    std::size_t raw_index = 0;
    std::streamoff offset = 0;
    double Ecm = std::numeric_limits<double>::quiet_NaN();
    int total_dim = 0;
    int proj_dim = 0;
    int success = 0;
};

struct RowResult {
    int row_id = -1;
    int cache_index = -1;
    double Ecm = std::numeric_limits<double>::quiet_NaN();
    double Lbyas = std::numeric_limits<double>::quiet_NaN();
    std::string irrep;
    int Nfull_gpu = 0;
    int Nproj_gpu = 0;
    int Nfull_cpu = 0;
    int Nproj_cpu = 0;
    comp det_gpu = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    comp det_cpu = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    double abs_diff = std::numeric_limits<double>::quiet_NaN();
    double rel_diff = std::numeric_limits<double>::quiet_NaN();
    std::string status = "OK";
};

struct BlockSummary {
    std::size_t raw_rows = 0;
    std::size_t unique_rows = 0;
    std::size_t selected_rows = 0;
    std::size_t duplicate_rows_removed = 0;
    std::size_t mismatch_rows = 0;
    std::size_t gpu_fail_rows = 0;
    std::size_t cpu_fail_rows = 0;
    double max_abs_diff = 0.0;
    double max_abs_diff_Ecm = std::numeric_limits<double>::quiet_NaN();
    double mean_abs_diff = std::numeric_limits<double>::quiet_NaN();
    double median_abs_diff = std::numeric_limits<double>::quiet_NaN();
    double max_rel_diff = 0.0;
    double max_rel_diff_Ecm = std::numeric_limits<double>::quiet_NaN();
    double mean_rel_diff = std::numeric_limits<double>::quiet_NaN();
    double rms_abs_diff = std::numeric_limits<double>::quiet_NaN();
    double avg_abs_det_gpu = std::numeric_limits<double>::quiet_NaN();
    double avg_abs_det_cpu = std::numeric_limits<double>::quiet_NaN();
    double avg_raw_det_gpu_real = std::numeric_limits<double>::quiet_NaN();
    double avg_raw_det_cpu_real = std::numeric_limits<double>::quiet_NaN();
    double max_abs_imag_gpu = std::numeric_limits<double>::quiet_NaN();
    double max_abs_imag_cpu = std::numeric_limits<double>::quiet_NaN();
};

static std::string trim(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    while(!s.empty() && is_space(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
    while(!s.empty() && is_space(static_cast<unsigned char>(s.back()))) s.pop_back();
    return s;
}

static std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for(char ch : s) {
        if(ch == ',') {
            if(!cur.empty()) out.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    if(!cur.empty()) out.push_back(trim(cur));
    return out;
}

static std::string ltag(double L) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(6) << L;
    std::string s = os.str();
    while(!s.empty() && s.back() == '0') s.pop_back();
    if(!s.empty() && s.back() == '.') s.pop_back();
    for(char& c : s) if(c == '.') c = 'p';
    return s;
}

static std::string block_prefix(double L, const std::string& irrep) {
    return std::string("L") + ltag(L) + "_" + irrep + "_oldscale_det_compare";
}

static void die_usage(const char* argv0) {
    std::cerr
        << "Usage: " << argv0 << " [options]\n"
        << "  --gpu-cache-dir PATH   root dir or cache file for coarseN=20000 GPU caches\n"
        << "  --outdir PATH          output directory\n"
        << "  --Lbyas 20,24          optional comma-separated L list\n"
        << "  --irrep 000_A1m,...    optional comma-separated irrep list\n"
        << "  --points N             selected rows per block (default 1000)\n"
        << "  --threads N            OMP threads (default 16)\n"
        << "  --k3df-zero            drop K3df and compare F3inv only\n"
        << "  --k3iso0 V --k3iso1 V --k3B V --k3E V\n";
}

static Args parse_args(int argc, char** argv) {
    Args a;
    for(int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto need_value = [&](const char* name) -> std::string {
            if(i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
            return argv[++i];
        };
        if(arg == "--help" || arg == "-h") {
            die_usage(argv[0]);
            std::exit(0);
        } else if(arg == "--gpu-cache-dir") {
            a.gpu_cache_dir = need_value("--gpu-cache-dir");
        } else if(arg == "--outdir") {
            a.outdir = need_value("--outdir");
        } else if(arg == "--plot-script") {
            a.plot_script = need_value("--plot-script");
        } else if(arg == "--Lbyas") {
            for(const auto& s : split_csv(need_value("--Lbyas"))) if(!s.empty()) a.Lbyas.push_back(std::stod(s));
        } else if(arg == "--irrep") {
            for(const auto& s : split_csv(need_value("--irrep"))) if(!s.empty()) a.irreps.push_back(s);
        } else if(arg == "--points") {
            a.points = std::max(1, std::stoi(need_value("--points")));
        } else if(arg == "--threads") {
            a.threads = std::max(1, std::stoi(need_value("--threads")));
        } else if(arg == "--xi") {
            a.xi = std::stod(need_value("--xi"));
        } else if(arg == "--k3df-zero") {
            a.k3df_zero = true;
        } else if(arg == "--k3iso0") {
            a.k3.K3iso0 = std::stod(need_value("--k3iso0"));
        } else if(arg == "--k3iso1") {
            a.k3.K3iso1 = std::stod(need_value("--k3iso1"));
        } else if(arg == "--k3B") {
            a.k3.K3B = std::stod(need_value("--k3B"));
        } else if(arg == "--k3E") {
            a.k3.K3E = std::stod(need_value("--k3E"));
        } else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }
    if(a.Lbyas.empty()) a.Lbyas = {20.0, 24.0};
    if(a.irreps.empty()) a.irreps = {"000_A1m", "100_A2", "110_A2", "111_A2", "200_A2"};
    if(a.k3df_zero) a.k3 = {0.0, 0.0, 0.0, 0.0};
    return a;
}

static std::string gpu_cache_filename(double L, const std::string& irrep) {
    std::ostringstream os;
    os << "v32zu_Lbyas" << ltag(L) << "_xi3p444_irrep" << irrep << "_coarse20000_F3inv_Vsel_gpu.bin";
    return os.str();
}

static std::optional<fs::path> resolve_gpu_cache(const fs::path& root, double L, const std::string& irrep) {
    if(root.empty()) return std::nullopt;
    if(fs::is_regular_file(root)) return root;
    const std::string want = gpu_cache_filename(L, irrep);
    fs::path direct = root / want;
    if(fs::is_regular_file(direct)) return direct;
    std::vector<fs::path> hits;
    if(fs::exists(root)) {
        for(auto it = fs::recursive_directory_iterator(root); it != fs::recursive_directory_iterator(); ++it) {
            if(!it->is_regular_file()) continue;
            if(it->path().filename() == want) hits.push_back(it->path());
        }
    }
    if(hits.empty()) return std::nullopt;
    std::sort(hits.begin(), hits.end());
    return hits.front();
}

static std::vector<GpuRowMeta> scan_gpu_cache(const fs::path& cache_path) {
    std::ifstream is(cache_path, std::ios::binary);
    if(!is) throw std::runtime_error("could not open GPU cache: " + cache_path.string());
    gpu_cache_reader_convention::require_gpu_cache_format_header(
        is,
        gpu_cache_reader_convention::GpuCacheFileFormat::CombinedRawF3invVselGpu,
        cache_path.string()
    );
    const std::uint64_t rec_magic = 0x5653325a4f524543ULL;
    std::vector<GpuRowMeta> rows;
    std::size_t idx = 0;
    while(true) {
        const std::streamoff offset = static_cast<std::streamoff>(is.tellg());
        if(!is.good()) break;
        std::uint64_t magic = 0;
        is.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        if(!is) break;
        if(magic != rec_magic) throw std::runtime_error("bad GPU cache record magic in " + cache_path.string());
        auto read_i32 = [&]() -> std::int32_t {
            std::int32_t x{};
            is.read(reinterpret_cast<char*>(&x), sizeof(x));
            if(!is) throw std::runtime_error("truncated GPU cache while scanning " + cache_path.string());
            return x;
        };
        auto read_f64 = [&]() -> double {
            double x{};
            is.read(reinterpret_cast<char*>(&x), sizeof(x));
            if(!is) throw std::runtime_error("truncated GPU cache while scanning " + cache_path.string());
            return x;
        };
        const std::int32_t grid_i = read_i32();
        const std::int32_t dim1 = read_i32();
        const std::int32_t dim2 = read_i32();
        const std::int32_t n = read_i32();
        const std::int32_t vdim = read_i32();
        const double Ecm = read_f64();
        const double En = read_f64();
        (void)grid_i;
        (void)dim1;
        (void)dim2;
        (void)En;
        const auto skip_matrix = [&](int rows, int cols) {
            const std::size_t bytes = std::size_t(rows) * std::size_t(cols) * 16u;
            is.seekg(static_cast<std::streamoff>(bytes), std::ios::cur);
            if(!is) throw std::runtime_error("failed to skip GPU matrix while scanning " + cache_path.string());
        };
        skip_matrix(n, n);
        skip_matrix(n, vdim);
        rows.push_back(GpuRowMeta{idx++, offset, Ecm, n, vdim, 1});
    }
    return rows;
}

static std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> load_gpu_row(const fs::path& cache_path,
                                                                   std::streamoff offset,
                                                                   v33g_runtime_k3basis::GpuComplexReadConvention conv = v33g_runtime_k3basis::GpuComplexReadConvention::RealImagSwappedV33dV32zu) {
    std::ifstream is(cache_path, std::ios::binary);
    if(!is) throw std::runtime_error("could not open GPU cache: " + cache_path.string());
    is.seekg(offset);
    if(!is) throw std::runtime_error("could not seek to GPU row offset in " + cache_path.string());
    auto read_i32 = [&]() -> std::int32_t {
        std::int32_t x{};
        is.read(reinterpret_cast<char*>(&x), sizeof(x));
        if(!is) throw std::runtime_error("truncated GPU cache row in " + cache_path.string());
        return x;
    };
    auto read_f64 = [&]() -> double {
        double x{};
        is.read(reinterpret_cast<char*>(&x), sizeof(x));
        if(!is) throw std::runtime_error("truncated GPU cache row in " + cache_path.string());
        return x;
    };
    std::uint64_t magic = 0;
    is.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    const std::uint64_t rec_magic = 0x5653325a4f524543ULL;
    if(magic != rec_magic) throw std::runtime_error("bad GPU cache record magic at seeked row in " + cache_path.string());
    const std::int32_t grid_i = read_i32();
    const std::int32_t dim1 = read_i32();
    const std::int32_t dim2 = read_i32();
    const std::int32_t n = read_i32();
    const std::int32_t vdim = read_i32();
    const double Ecm = read_f64();
    const double En = read_f64();
    (void)grid_i;
    (void)dim1;
    (void)dim2;
    (void)Ecm;
    (void)En;
    if(conv != v33g_runtime_k3basis::GpuComplexReadConvention::RealImagSwappedV33dV32zu) {
        throw std::runtime_error("unsupported GPU cache complex read convention in load_gpu_row");
    }
    Eigen::MatrixXcd F3inv = v33g_runtime_k3basis::read_gpu_cache_matrix_col_major(is, n, n, conv);
    Eigen::MatrixXcd Vsel = v33g_runtime_k3basis::read_gpu_cache_matrix_col_major(is, n, vdim, conv);
    return {std::move(F3inv), std::move(Vsel)};
}

static PhysicsParams make_physics(double L, double xi, int threads) {
    PhysicsParams par;
    par.atmpi = 0.06906;
    par.atmK = 0.09698;
    par.eta_1 = 1.0;
    par.eta_2 = 0.5;
    par.alpha = 0.5;
    par.epsilon_h = 0.0;
    par.max_shell_num = 20.0;
    par.tolerance = 1.0e-12;
    par.xi = xi;
    par.Lbyas = L;
    par.parity = -1;
    par.eig_tol = 0.05;
    par.norm_tol = 1.0e-12;
    par.proj_tol = 1.0e-10;
    par.omp_threads = threads;
    par.waves_vec_1 = {0, 1};
    par.waves_vec_2 = {0};
    par.scatter_params_1.assign(4, std::vector<comp>(3, comp(0.0, 0.0)));
    par.scatter_params_2.assign(4, std::vector<comp>(3, comp(0.0, 0.0)));
    par.scatter_params_1[0][0] = comp(4.04, 0.0);
    par.scatter_params_1[1][0] = comp(-43.2, 0.0);
    par.scatter_params_2[0][0] = comp(4.12, 0.0);
    return par;
}

static FitSettings make_settings(double L, double xi, int threads, const std::string& irrep) {
    FitSettings s;
    s.Lval = L;
    s.xival = xi;
    s.scan_E0 = 0.2631;
    s.scan_E1 = 0.36;
    s.coarseN = 20000;
    s.omp_threads = threads;
    s.debug = 'n';
    s.atmpi = 0.06906;
    s.atmK = 0.09698;
    s.eta_1 = 1.0;
    s.eta_2 = 0.5;
    s.alpha = 0.5;
    s.epsilon_h = 0.0;
    s.max_shell_num = 20.0;
    s.tolerance = 1.0e-12;
    s.parity = -1;
    s.eig_tol = 0.05;
    s.norm_tol = 1.0e-12;
    s.proj_tol = 1.0e-10;
    s.waves_vec_1 = {0, 1};
    s.waves_vec_2 = {0};
    s.scatter_params_1[0][0] = 4.04;
    s.scatter_params_1[1][0] = -43.2;
    s.scatter_params_2[0][0] = 4.12;
    s.list_of_mom = {irrep};
    return s;
}

static comp det_complex(const Eigen::MatrixXcd& M) {
    if(M.rows() <= 0 || M.rows() != M.cols()) {
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu(M);
    return lu.determinant();
}

static void write_csv_header(std::ofstream& os) {
    os << "row_id,cache_index,Ecm,Lbyas,irrep,Nfull_gpu,Nproj_gpu,Nfull_cpu,Nproj_cpu,"
       << "det_gpu_real,det_gpu_imag,det_cpu_real,det_cpu_imag,abs_diff,rel_diff,status\n";
}

static BlockSummary summarize(const std::vector<RowResult>& rows) {
    BlockSummary s;
    std::vector<double> abs_vals;
    abs_vals.reserve(rows.size());
    double sum_abs = 0.0;
    double sum_rel = 0.0;
    double sum_abs2 = 0.0;
    double sum_abs_det_gpu = 0.0;
    double sum_abs_det_cpu = 0.0;
    double sum_raw_det_gpu_real = 0.0;
    double sum_raw_det_cpu_real = 0.0;
    double max_abs_imag_gpu = 0.0;
    double max_abs_imag_cpu = 0.0;
    std::size_t n_ok = 0;
    for(const auto& r : rows) {
        if(r.status == "DIM_MISMATCH") ++s.mismatch_rows;
        if(r.status == "GPU_FAIL" || r.status == "NUMERIC_FAIL") ++s.gpu_fail_rows;
        if(r.status == "CPU_FAIL") ++s.cpu_fail_rows;
        if(r.status != "OK") continue;
        ++n_ok;
        abs_vals.push_back(r.abs_diff);
        sum_abs += r.abs_diff;
        sum_rel += r.rel_diff;
        sum_abs2 += r.abs_diff * r.abs_diff;
        sum_abs_det_gpu += std::abs(r.det_gpu);
        sum_abs_det_cpu += std::abs(r.det_cpu);
        sum_raw_det_gpu_real += r.det_gpu.real();
        sum_raw_det_cpu_real += r.det_cpu.real();
        max_abs_imag_gpu = std::max(max_abs_imag_gpu, std::abs(r.det_gpu.imag()));
        max_abs_imag_cpu = std::max(max_abs_imag_cpu, std::abs(r.det_cpu.imag()));
        if(r.abs_diff > s.max_abs_diff) {
            s.max_abs_diff = r.abs_diff;
            s.max_abs_diff_Ecm = r.Ecm;
        }
        if(r.rel_diff > s.max_rel_diff) {
            s.max_rel_diff = r.rel_diff;
            s.max_rel_diff_Ecm = r.Ecm;
        }
    }
    s.mean_abs_diff = n_ok ? sum_abs / double(n_ok) : std::numeric_limits<double>::quiet_NaN();
    s.mean_rel_diff = n_ok ? sum_rel / double(n_ok) : std::numeric_limits<double>::quiet_NaN();
    s.rms_abs_diff = n_ok ? std::sqrt(sum_abs2 / double(n_ok)) : std::numeric_limits<double>::quiet_NaN();
    s.avg_abs_det_gpu = n_ok ? sum_abs_det_gpu / double(n_ok) : std::numeric_limits<double>::quiet_NaN();
    s.avg_abs_det_cpu = n_ok ? sum_abs_det_cpu / double(n_ok) : std::numeric_limits<double>::quiet_NaN();
    s.avg_raw_det_gpu_real = n_ok ? sum_raw_det_gpu_real / double(n_ok) : std::numeric_limits<double>::quiet_NaN();
    s.avg_raw_det_cpu_real = n_ok ? sum_raw_det_cpu_real / double(n_ok) : std::numeric_limits<double>::quiet_NaN();
    s.max_abs_imag_gpu = max_abs_imag_gpu;
    s.max_abs_imag_cpu = max_abs_imag_cpu;
    if(!abs_vals.empty()) {
        auto mid = abs_vals.begin() + abs_vals.size() / 2;
        std::nth_element(abs_vals.begin(), mid, abs_vals.end());
        s.median_abs_diff = *mid;
    }
    return s;
}

static void print_block_summary(std::ostream& os, const BlockSummary& s) {
    os << std::scientific << std::setprecision(17);
    os << "[summary] mismatch_count = " << s.mismatch_rows << "\n";
    os << "[summary] gpu_fail_count = " << s.gpu_fail_rows << "\n";
    os << "[summary] cpu_fail_count = " << s.cpu_fail_rows << "\n";
    os << "[summary] max_abs_diff = " << s.max_abs_diff << " at Ecm=" << s.max_abs_diff_Ecm << "\n";
    os << "[summary] mean_abs_diff = " << s.mean_abs_diff << "\n";
    os << "[summary] median_abs_diff = " << s.median_abs_diff << "\n";
    os << "[summary] max_rel_diff = " << s.max_rel_diff << " at Ecm=" << s.max_rel_diff_Ecm << "\n";
    os << "[summary] mean_rel_diff = " << s.mean_rel_diff << "\n";
    os << "[summary] rms_abs_diff = " << s.rms_abs_diff << "\n";
    os << "[summary] avg_abs_det_gpu = " << s.avg_abs_det_gpu << "\n";
    os << "[summary] avg_abs_det_cpu = " << s.avg_abs_det_cpu << "\n";
    os << "[summary] avg_raw_det_gpu_real = " << s.avg_raw_det_gpu_real << "\n";
    os << "[summary] avg_raw_det_cpu_real = " << s.avg_raw_det_cpu_real << "\n";
    os << "[summary] max_abs_imag_gpu = " << s.max_abs_imag_gpu << "\n";
    os << "[summary] max_abs_imag_cpu = " << s.max_abs_imag_cpu << "\n";
}

static std::string build_param_block(double L, double xi, const K3dfParameters& k3, bool zero_mode) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(17);
    os << "[params] xi=" << xi << " Lbyas=" << L << " Lxi=" << (L * xi) << "\n";
    os << "[params] atmK=0.09698 atmpi=0.06906 eta_1=1 eta_2=0.5 alpha=0.5 epsilon_h=0 max_shell_num=20 tolerance=1e-12 parity=-1 eig_tol=0.05 norm_tol=1e-12 proj_tol=1e-10\n";
    os << "[params] waves_vec_1={0,1} waves_vec_2={0}\n";
    os << "[params] scatter1_00=4.04 scatter1_10=-43.2 scatter2_00=4.12\n";
    os << "[params] k3df_zero=" << (zero_mode ? "true" : "false")
       << " K3iso0=" << k3.K3iso0
       << " K3iso1=" << k3.K3iso1
       << " K3B=" << k3.K3B
       << " K3E=" << k3.K3E << "\n";
    return os.str();
}

static int run_plot_script(const fs::path& script, const fs::path& csv, const fs::path& outdir, const std::string& prefix) {
    std::ostringstream cmd;
    cmd << "python3 " << std::quoted(script.string())
        << " --csv " << std::quoted(csv.string())
        << " --outdir " << std::quoted(outdir.string())
        << " --prefix " << std::quoted(prefix);
    const int rc = std::system(cmd.str().c_str());
    return rc;
}

} // namespace

#ifndef DISABLE_OLDSCALE_COMPARE_MAIN
int main(int argc, char** argv) {
    try {
        const Args args = parse_args(argc, argv);
        fs::create_directories(args.outdir);
        if(args.threads > 0) {
#ifdef _OPENMP
            omp_set_dynamic(0);
            omp_set_num_threads(args.threads);
            omp_set_max_active_levels(1);
#endif
            Eigen::setNbThreads(1);
        }

        std::cout << std::scientific << std::setprecision(17);
        std::cout << "[config] gpu_cache_dir=" << args.gpu_cache_dir << "\n";
        std::cout << "[config] outdir=" << args.outdir << "\n";
        std::cout << "[config] plot_script=" << args.plot_script << "\n";
        std::cout << "[config] points=" << args.points << " threads=" << args.threads << "\n";
        std::cout << "[config] Lbyas_list=";
        for(std::size_t i = 0; i < args.Lbyas.size(); ++i) {
            if(i) std::cout << ",";
            std::cout << args.Lbyas[i];
        }
        std::cout << "\n[config] irreps_list=";
        for(std::size_t i = 0; i < args.irreps.size(); ++i) {
            if(i) std::cout << ",";
            std::cout << args.irreps[i];
        }
        std::cout << "\n";

        std::vector<std::string> block_reports;
        for(double L : args.Lbyas) {
            for(const auto& irrep : args.irreps) {
                const auto cache_path_opt = resolve_gpu_cache(args.gpu_cache_dir, L, irrep);
                const fs::path cache_path = cache_path_opt ? *cache_path_opt : fs::path();
                if(cache_path.empty() || !fs::exists(cache_path)) {
                    std::cout << "[warn] missing GPU cache for Lbyas=" << L << " irrep=" << irrep << " under " << args.gpu_cache_dir << "\n";
                    continue;
                }
                const fs::path meta_path = cache_path.string() + ".meta.json";
                if(!fs::exists(meta_path)) {
                    std::cout << "[warn] missing meta sidecar, skipping block: " << cache_path << "\n";
                    continue;
                }
                const auto meta = v33g_runtime_k3basis::read_runtime_meta_kv(cache_path.string());
                const std::string version = v33g_runtime_k3basis::meta_string(meta, "version", "");
                const std::string cache_kind = v33g_runtime_k3basis::meta_string(meta, "cache_kind", "");
                const double meta_L = v33g_runtime_k3basis::meta_double(meta, "Lbyas", -1.0);
                const double meta_xi = v33g_runtime_k3basis::meta_double(meta, "xi", -1.0);
                const std::string meta_irrep = v33g_runtime_k3basis::meta_string(meta, "irrep_label", v33g_runtime_k3basis::meta_string(meta, "irrep", ""));
                const std::string meta_irrep_tag = v33g_runtime_k3basis::meta_string(meta, "irrep_tag", meta_irrep);
                const int meta_coarseN = v33g_runtime_k3basis::meta_int(meta, "coarseN", -1);
                const int meta_rows = v33g_runtime_k3basis::meta_int(meta, "rows_written", v33g_runtime_k3basis::meta_int(meta, "rows", -1));
                const bool version_ok = (version == "v33d_gpu_component_cachegen_all_ingredients" || version == "v33g_runtime_k3basis_cache");
                const bool kind_ok = (cache_kind == "F3inv_plus_Vsel_gpu_full" || cache_kind == "v33g_runtime_k3basis_cache" || cache_kind.empty());
                if(!version_ok || !kind_ok || std::abs(meta_L - L) > 1e-12 || std::abs(meta_xi - args.xi) > 1e-12 || (meta_irrep != irrep && meta_irrep_tag != irrep) || meta_coarseN != 20000 || meta_rows != 20000) {
                    std::cout << "[warn] cache metadata mismatch, skipping block: " << cache_path
                              << " version=" << version << " cache_kind=" << cache_kind
                              << " L=" << meta_L << " xi=" << meta_xi << " irrep_label=" << meta_irrep
                              << " irrep_tag=" << meta_irrep_tag
                              << " coarseN=" << meta_coarseN << " rows=" << meta_rows << "\n";
                    continue;
                }

                const std::vector<GpuRowMeta> raw_rows = scan_gpu_cache(cache_path);
                if(raw_rows.empty()) {
                    std::cout << "[warn] no GPU rows found in " << cache_path << "\n";
                    continue;
                }
                std::vector<GpuRowMeta> rows = raw_rows;
                std::sort(rows.begin(), rows.end(), [](const GpuRowMeta& a, const GpuRowMeta& b) {
                    if(std::abs(a.Ecm - b.Ecm) > 1e-13) return a.Ecm < b.Ecm;
                    return a.raw_index < b.raw_index;
                });
                std::vector<GpuRowMeta> unique_rows;
                unique_rows.reserve(rows.size());
                std::size_t duplicate_rows_removed = 0;
                for(const auto& row : rows) {
                    if(unique_rows.empty() || std::abs(row.Ecm - unique_rows.back().Ecm) > 1e-13) {
                        unique_rows.push_back(row);
                    } else {
                        ++duplicate_rows_removed;
                    }
                }
                const std::size_t available = unique_rows.size();
                std::size_t n_select = args.points;
                if(available < n_select) {
                    std::cout << "[warn] requested points=" << n_select << " but only " << available << " unique rows exist; using all available rows\n";
                    n_select = available;
                }
                std::vector<std::size_t> pick_indices;
                pick_indices.reserve(n_select);
                if(n_select == 1) {
                    pick_indices.push_back(0);
                } else {
                    for(std::size_t j = 0; j < n_select; ++j) {
                        const std::size_t idx = static_cast<std::size_t>(std::llround(double(j) * double(available - 1) / double(n_select - 1)));
                        pick_indices.push_back(std::min(idx, available - 1));
                    }
                }
                std::vector<std::pair<std::size_t, GpuRowMeta>> selected;
                selected.reserve(n_select);
                for(std::size_t idx : pick_indices) selected.push_back({idx, unique_rows[idx]});

                const double const_norm_scale = std::pow(L * args.xi, 6.0);
                const std::string prefix = block_prefix(L, irrep);
                const fs::path csv_path = args.outdir / (prefix + ".csv");
                const fs::path summary_path = args.outdir / (prefix + "_summary.txt");
                std::ofstream csv(csv_path);
                if(!csv) throw std::runtime_error("could not open CSV for writing: " + csv_path.string());
                csv << std::scientific << std::setprecision(17);
                write_csv_header(csv);

                std::cout << "[block] Lbyas=" << L << " irrep=" << irrep
                          << " cache_rows=" << raw_rows.size()
                          << " selected_points=" << selected.size() << "\n";
                std::cout << "[block] unique_rows=" << unique_rows.size()
                          << " duplicate_rows_removed=" << duplicate_rows_removed << "\n";
                std::cout << "[scaling] old const_norm_scale=(Lbyas*xi)^6 = " << const_norm_scale << "\n";
                std::cout << "[det] det = det( V^dagger * ((F3inv+K3df)/const_norm_scale) * V )\n";
                std::cout << build_param_block(L, args.xi, args.k3, args.k3df_zero);

                const PhysicsParams par = make_physics(L, args.xi, args.threads);
                const FitSettings settings = make_settings(L, args.xi, args.threads, irrep);
                std::vector<RowResult> results(selected.size());

                #pragma omp parallel for schedule(dynamic, 1)
                for(int sidx = 0; sidx < static_cast<int>(selected.size()); ++sidx) {
                    const auto selected_row = selected[std::size_t(sidx)];
                    const std::size_t cache_index = selected_row.first;
                    const GpuRowMeta row = selected_row.second;
                    RowResult rr;
                    rr.row_id = sidx;
                    rr.cache_index = static_cast<int>(cache_index);
                    rr.Ecm = row.Ecm;
                    rr.Lbyas = L;
                    rr.irrep = irrep;
                    rr.Nfull_gpu = row.total_dim;
                    rr.Nproj_gpu = row.proj_dim;
                    rr.status = "OK";

                    try {
                        const auto [F3inv_gpu, Vsel_gpu] = load_gpu_row(cache_path, row.offset);
                        rr.Nfull_gpu = F3inv_gpu.rows();
                        rr.Nproj_gpu = Vsel_gpu.cols();
                        ProjectedQCCacheEntry cpu_entry = build_cache_entry(static_cast<int>(cache_index), row.Ecm, parse_label(irrep), settings, par, 'n');
                        rr.Nfull_cpu = cpu_entry.total_dim;
                        rr.Nproj_cpu = cpu_entry.proj_dim;
                        if(!cpu_entry.success) {
                            rr.status = "CPU_FAIL";
                            rr.det_gpu = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
                            rr.det_cpu = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
                            rr.abs_diff = std::numeric_limits<double>::quiet_NaN();
                            rr.rel_diff = std::numeric_limits<double>::quiet_NaN();
                            results[std::size_t(sidx)] = rr;
                            continue;
                        }
                        if(rr.Nfull_gpu != rr.Nfull_cpu || rr.Nproj_gpu != rr.Nproj_cpu) {
                            rr.status = "DIM_MISMATCH";
                        }

                        if(rr.status == "OK") {
                            Eigen::MatrixXcd K3_full(cpu_entry.total_dim, cpu_entry.total_dim);
                            std::vector<comp> Kiso = {comp(args.k3.K3iso0, 0.0), comp(args.k3.K3iso1, 0.0)};
                            k3_2plus1::K3mat_2plus1(
                                K3_full,
                                cpu_entry.En_c,
                                cpu_entry.plm_config,
                                cpu_entry.klm_config,
                                cpu_entry.total_P,
                                par.atmK,
                                par.atmpi,
                                Kiso,
                                comp(args.k3.K3B, 0.0),
                                comp(args.k3.K3E, 0.0),
                                'n'
                            );
                            Eigen::MatrixXcd qc_gpu = F3inv_gpu + K3_full;
                            Eigen::MatrixXcd qc_cpu = cpu_entry.F3inv_full + K3_full;
                            qc_gpu = qc_gpu / const_norm_scale;
                            qc_cpu = qc_cpu / const_norm_scale;
                            Eigen::MatrixXcd proj_gpu = Vsel_gpu.adjoint() * qc_gpu * Vsel_gpu;
                            Eigen::MatrixXcd proj_cpu = cpu_entry.Vsel.adjoint() * qc_cpu * cpu_entry.Vsel;
                            rr.det_gpu = det_complex(proj_gpu);
                            rr.det_cpu = det_complex(proj_cpu);
                            rr.abs_diff = std::abs(rr.det_gpu - rr.det_cpu);
                            const double tiny = 1.0e-300;
                            rr.rel_diff = rr.abs_diff / std::max({std::abs(rr.det_gpu), std::abs(rr.det_cpu), tiny});
                            if(!std::isfinite(rr.abs_diff) || !std::isfinite(rr.rel_diff)) rr.status = "NUMERIC_FAIL";
                        }
                    } catch(const std::exception& ex) {
                        rr.status = "GPU_FAIL";
                        rr.det_gpu = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
                        rr.det_cpu = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
                        rr.abs_diff = std::numeric_limits<double>::quiet_NaN();
                        rr.rel_diff = std::numeric_limits<double>::quiet_NaN();
                        rr.Nfull_cpu = 0;
                        rr.Nproj_cpu = 0;
                        std::cerr << "[row-error] " << ex.what() << "\n";
                    }
                    results[std::size_t(sidx)] = rr;
                }

                std::sort(results.begin(), results.end(), [](const RowResult& a, const RowResult& b) {
                    return a.row_id < b.row_id;
                });
                for(const auto& rr : results) {
                    csv << rr.row_id << "," << rr.cache_index << "," << rr.Ecm << "," << rr.Lbyas << "," << rr.irrep << ","
                        << rr.Nfull_gpu << "," << rr.Nproj_gpu << "," << rr.Nfull_cpu << "," << rr.Nproj_cpu << ","
                        << rr.det_gpu.real() << "," << rr.det_gpu.imag() << ","
                        << rr.det_cpu.real() << "," << rr.det_cpu.imag() << ","
                        << rr.abs_diff << "," << rr.rel_diff << "," << rr.status << "\n";
                    std::cout << std::setw(5) << rr.row_id << "  "
                              << std::setprecision(17) << rr.Ecm << "  "
                              << rr.Nfull_gpu << "  " << rr.Nproj_gpu << "  "
                              << rr.Nfull_cpu << "  " << rr.Nproj_cpu << "  "
                              << std::showpos << std::scientific
                              << rr.det_gpu.real() << "  " << rr.det_gpu.imag() << "  "
                              << rr.det_cpu.real() << "  " << rr.det_cpu.imag() << "  "
                              << std::noshowpos
                              << rr.abs_diff << "  " << rr.rel_diff << "  "
                              << "STATUS=" << rr.status << "\n";
                }

                const BlockSummary summary = summarize(results);
                std::ostringstream sumtxt;
                sumtxt << std::scientific << std::setprecision(17);
                sumtxt << "[block] Lbyas=" << L << " irrep=" << irrep << "\n";
                sumtxt << "[block] cache_file=" << cache_path << "\n";
                sumtxt << "[block] raw_rows=" << raw_rows.size() << " unique_rows=" << unique_rows.size()
                       << " selected_rows=" << selected.size() << " duplicate_rows_removed=" << duplicate_rows_removed << "\n";
                sumtxt << "[scaling] old const_norm_scale=(Lbyas*xi)^6 = " << const_norm_scale << "\n";
                sumtxt << build_param_block(L, args.xi, args.k3, args.k3df_zero);
                print_block_summary(sumtxt, summary);
                if(std::max(summary.max_abs_imag_gpu, summary.max_abs_imag_cpu) > 1.0e-10) {
                    sumtxt << "[warning] determinant imaginary part is non-negligible; overlay plot shows real parts only\n";
                }
                if(summary.mismatch_rows > 0) {
                    sumtxt << "[warning] " << summary.mismatch_rows << " rows had dimension mismatch and were not compared\n";
                }
                csv.close();
                {
                    std::ofstream out(summary_path);
                    out << sumtxt.str();
                }
                std::cout << sumtxt.str();

                const int plot_rc = run_plot_script(args.plot_script, csv_path, args.outdir, prefix);
                if(plot_rc != 0) {
                    throw std::runtime_error("plot script failed for " + csv_path.string());
                }
                block_reports.push_back(summary_path.string());
            }
        }
        if(block_reports.empty()) {
            std::cout << "[done] no blocks were processed\n";
            return 0;
        }
        std::cout << "[done] processed " << block_reports.size() << " block(s)\n";
        return 0;
    } catch(const std::exception& ex) {
        std::cerr << "[fatal] " << ex.what() << "\n";
        return 1;
    }
}
#endif
