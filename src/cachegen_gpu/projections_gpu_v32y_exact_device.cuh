#ifndef PROJECTIONS_GPU_V32Y_EXACT_DEVICE_CUH
#define PROJECTIONS_GPU_V32Y_EXACT_DEVICE_CUH

// v32zt exact-v32y GPU Vsel projector.
// This is the GPU-only production projector used by the cache generator.
// It ports the v32y stream-cache Vsel path validated by v32zs:
//   config_maker_4_momentum_first basis ordering
//   fvproj_v30j::projector_for_convention convention index=0
//   use_inverse_momentum=false, use_inverse_D=false, transpose_D=false,
//   order_WS=false, parity_mode=1, global parity=-1
//   v30m Cartesian ell=1 action: Y_11~x, Y_10~z, Y_1,-1~y
//   Vsel = eigenvectors of P_I with eigenvalues near 1, eig_tol=0.05.

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cuComplex.h>
#include <vector>
#include <array>
#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>

namespace v32zt_vsel_gpu {

__device__ inline int det_sign_perm_dev(const int* R)
{
    int p0 = abs(R[0]) - 1;
    int p1 = abs(R[1]) - 1;
    int p2 = abs(R[2]) - 1;
    int inv = (p0 > p1) + (p0 > p2) + (p1 > p2);
    int psign = (inv % 2 == 0) ? 1 : -1;
    int s = (R[0] > 0 ? 1 : -1) * (R[1] > 0 ? 1 : -1) * (R[2] > 0 ? 1 : -1);
    return psign * s;
}

__device__ inline int cart_axis_from_m_dev(int m)
{
    if (m ==  1) return 0; // x
    if (m ==  0) return 2; // z
    if (m == -1) return 1; // y
    return -1;
}

__device__ inline double D_real_v32y_v30m_dev(int ell, int m_tgt, int m_src, const int* R)
{
    if (ell == 0) return (m_tgt == 0 && m_src == 0) ? 1.0 : 0.0;
    if (ell == 1) {
        int at = cart_axis_from_m_dev(m_tgt);
        int as = cart_axis_from_m_dev(m_src);
        if (at < 0 || as < 0) return 0.0;
        int src_axis_for_target = abs(R[at]) - 1;
        int sign = (R[at] > 0) ? 1 : -1;
        return (src_axis_for_target == as) ? double(sign) : 0.0;
    }
    // Production default is waves_vec_1={0,1}, waves_vec_2={0}.
    // Returning 0 for ell>1 prevents silent unsupported behavior.
    return 0.0;
}

__device__ inline int transform_component_dev(int row, int sx, int sy, int sz, const int* R)
{
    int r = R[row];
    int idx = abs(r) - 1;
    int sign = (r > 0) ? 1 : -1;
    int v = (idx == 0 ? sx : (idx == 1 ? sy : sz));
    return sign * v;
}

__global__ void build_projector_v32y_exact_kernel(
    int N,
    int nrot,
    int dI,
    int parity,
    const int* nx,
    const int* ny,
    const int* nz,
    const int* ell,
    const int* mm,
    const int* flavor,
    const int* Rflat,
    const double* chivals,
    double* P)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = N * N;
    if (idx >= total) return;

    int target = idx % N;
    int source = idx / N; // column-major P[target + source*N]
    double val = 0.0;

    if (flavor[target] == flavor[source] && ell[target] == ell[source]) {
        int es = ell[source];
        int mt = mm[target];
        int ms = mm[source];
        int sx = nx[source];
        int sy = ny[source];
        int sz = nz[source];

        for (int r = 0; r < nrot; ++r) {
            const int* R = Rflat + 3*r;
            int tx = transform_component_dev(0, sx, sy, sz, R);
            int ty = transform_component_dev(1, sx, sy, sz, R);
            int tz = transform_component_dev(2, sx, sy, sz, R);
            if (tx != nx[target] || ty != ny[target] || tz != nz[target]) continue;

            double pf = 1.0;
            if (parity == -1 && det_sign_perm_dev(R) < 0) pf = -1.0; // v32y parity_mode=1 convention
            double d = D_real_v32y_v30m_dev(es, mt, ms, R);
            val += chivals[r] * pf * d;
        }
    }

    if (nrot > 0) val *= (double(dI) / double(nrot));
    if (fabs(val) < 1.0e-15) val = 0.0;
    P[target + source*N] = val;
}

__global__ void symmetrize_real_matrix_v32y_kernel(int N, double* P)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = N * N;
    if (idx >= total) return;
    int i = idx % N;
    int j = idx / N;
    if (i < j) {
        double a = P[i + j*N];
        double b = P[j + i*N];
        double c = 0.5 * (a + b);
        P[i + j*N] = c;
        P[j + i*N] = c;
    }
}

__global__ void gather_v32y_evecs_to_complex_kernel(
    cuDoubleComplex* Vsel,
    const double* evecs,
    const int* selected,
    int N,
    int vdim,
    double chop_tol)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = N * vdim;
    if (tid >= total) return;
    int row = tid % N;
    int col = tid / N;
    int src = selected[col];
    double v = evecs[row + src * N];
    if (fabs(v) < chop_tol) v = 0.0;
    Vsel[tid] = make_cuDoubleComplex(v, 0.0);
}

inline void build_basis_arrays_exact_v32y(
    const std::vector<std::vector<int>>& np_config,
    const std::vector<std::vector<int>>& nk_config,
    std::vector<int>& nx,
    std::vector<int>& ny,
    std::vector<int>& nz,
    std::vector<int>& ell,
    std::vector<int>& m,
    std::vector<int>& flavor)
{
    if (np_config.size() < 5 || nk_config.size() < 5)
        throw std::runtime_error("v32y exact Vsel: n_config rows < 5");
    const int n1 = static_cast<int>(np_config[0].size());
    const int n2 = static_cast<int>(nk_config[0].size());
    const int N = n1 + n2;
    nx.resize(N); ny.resize(N); nz.resize(N); ell.resize(N); m.resize(N); flavor.resize(N);
    for (int i=0; i<n1; ++i) {
        nx[i] = np_config[0][i]; ny[i] = np_config[1][i]; nz[i] = np_config[2][i];
        ell[i] = np_config[3][i]; m[i] = np_config[4][i]; flavor[i] = 1;
    }
    for (int j=0; j<n2; ++j) {
        const int i = n1 + j;
        nx[i] = nk_config[0][j]; ny[i] = nk_config[1][j]; nz[i] = nk_config[2][j];
        ell[i] = nk_config[3][j]; m[i] = nk_config[4][j]; flavor[i] = 2;
    }
}

inline void copy_int_vector_to_device_exact_v32y(const std::vector<int>& h, int*& d)
{
    V11_CUDA_CHECK(cudaMalloc(&d, sizeof(int) * h.size()));
    V11_CUDA_CHECK(cudaMemcpy(d, h.data(), sizeof(int) * h.size(), cudaMemcpyHostToDevice));
}

inline void build_projection_Vsel_device_exact_v32y(
    cuDoubleComplex*& dVsel_out,
    int& vdim_out,
    int N,
    int size1,
    const std::vector<std::vector<int>>& np_config,
    const std::vector<std::vector<int>>& nk_config,
    const std::string& irrep,
    int nnP0,
    int nnP1,
    int nnP2,
    int parity,
    double eig_tol,
    double chop_tol,
    int threads,
    char debug)
{
    dVsel_out = nullptr;
    vdim_out = 0;
    if (N <= 0) return;

    std::vector<int> h_nx, h_ny, h_nz, h_ell, h_m, h_flavor;
    build_basis_arrays_exact_v32y(np_config, nk_config, h_nx, h_ny, h_nz, h_ell, h_m, h_flavor);
    if (static_cast<int>(h_nx.size()) != N)
        throw std::runtime_error("v32y exact GPU projection: basis size mismatch");
    if (size1 != static_cast<int>(np_config[0].size()))
        throw std::runtime_error("v32y exact GPU projection: size1 mismatch");

    auto LG = pgpu::little_group_host(nnP0, nnP1, nnP2);
    const int nrot = static_cast<int>(LG.size());
    if (nrot <= 0) throw std::runtime_error("v32y exact GPU projection: empty little group");
    const int dI = pgpu::irrep_dim_host(irrep);

    std::vector<int> hR;
    hR.reserve(static_cast<std::size_t>(3*nrot));
    std::vector<double> hchi;
    hchi.reserve(static_cast<std::size_t>(nrot));
    for (const auto& R : LG) {
        hR.push_back(R[0]); hR.push_back(R[1]); hR.push_back(R[2]);
        hchi.push_back(static_cast<double>(pgpu::chi_host(R, irrep, nnP0, nnP1, nnP2)));
    }

    int *d_nx=nullptr, *d_ny=nullptr, *d_nz=nullptr, *d_ell=nullptr, *d_m=nullptr, *d_flavor=nullptr, *d_R=nullptr, *d_sel=nullptr, *d_info=nullptr;
    double *d_chi=nullptr, *d_P=nullptr, *d_W=nullptr, *d_work=nullptr;
    cusolverDnHandle_t solver=nullptr;

    try {
        copy_int_vector_to_device_exact_v32y(h_nx, d_nx);
        copy_int_vector_to_device_exact_v32y(h_ny, d_ny);
        copy_int_vector_to_device_exact_v32y(h_nz, d_nz);
        copy_int_vector_to_device_exact_v32y(h_ell, d_ell);
        copy_int_vector_to_device_exact_v32y(h_m, d_m);
        copy_int_vector_to_device_exact_v32y(h_flavor, d_flavor);
        V11_CUDA_CHECK(cudaMalloc(&d_R, sizeof(int) * hR.size()));
        V11_CUDA_CHECK(cudaMemcpy(d_R, hR.data(), sizeof(int)*hR.size(), cudaMemcpyHostToDevice));
        V11_CUDA_CHECK(cudaMalloc(&d_chi, sizeof(double) * hchi.size()));
        V11_CUDA_CHECK(cudaMemcpy(d_chi, hchi.data(), sizeof(double)*hchi.size(), cudaMemcpyHostToDevice));
        V11_CUDA_CHECK(cudaMalloc(&d_P, sizeof(double) * static_cast<std::size_t>(N) * static_cast<std::size_t>(N)));

        const int blocks = (N*N + threads - 1) / threads;
        build_projector_v32y_exact_kernel<<<blocks, threads>>>(N, nrot, dI, parity,
            d_nx, d_ny, d_nz, d_ell, d_m, d_flavor, d_R, d_chi, d_P);
        V11_CUDA_CHECK(cudaGetLastError());
        symmetrize_real_matrix_v32y_kernel<<<blocks, threads>>>(N, d_P);
        V11_CUDA_CHECK(cudaGetLastError());
        V11_CUDA_CHECK(cudaDeviceSynchronize());

        V11_CUSOLVER_CHECK(cusolverDnCreate(&solver));
        V11_CUDA_CHECK(cudaMalloc(&d_W, sizeof(double)*N));
        V11_CUDA_CHECK(cudaMalloc(&d_info, sizeof(int)));
        int lwork = 0;
        V11_CUSOLVER_CHECK(cusolverDnDsyevd_bufferSize(solver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
                                                       N, d_P, N, d_W, &lwork));
        V11_CUDA_CHECK(cudaMalloc(&d_work, sizeof(double)*static_cast<std::size_t>(lwork)));
        V11_CUSOLVER_CHECK(cusolverDnDsyevd(solver, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,
                                            N, d_P, N, d_W, d_work, lwork, d_info));
        V11_CUDA_CHECK(cudaDeviceSynchronize());
        int info = 0;
        V11_CUDA_CHECK(cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
        if (info != 0) throw std::runtime_error("v32y exact GPU projection Dsyevd failed info=" + std::to_string(info));

        std::vector<double> evals(static_cast<std::size_t>(N));
        V11_CUDA_CHECK(cudaMemcpy(evals.data(), d_W, sizeof(double)*evals.size(), cudaMemcpyDeviceToHost));
        std::vector<int> selected;
        selected.reserve(N);
        for (int i=0; i<N; ++i) {
            const double ev = evals[static_cast<std::size_t>(i)];
            if (ev >= 1.0 - eig_tol && ev <= 1.0 + eig_tol) selected.push_back(i);
        }
        vdim_out = static_cast<int>(selected.size());
        if (vdim_out <= 0) {
            std::ostringstream os;
            os << "v32y exact GPU projection produced Vsel with zero columns"
               << " N=" << N << " size1=" << size1 << " irrep=" << irrep
               << " P=[" << nnP0 << "," << nnP1 << "," << nnP2 << "]"
               << " parity=" << parity << " eig_tol=" << eig_tol
               << " eval_min=" << (N ? evals.front() : 0.0)
               << " eval_max=" << (N ? evals.back() : 0.0);
            throw std::runtime_error(os.str());
        }

        V11_CUDA_CHECK(cudaMalloc(&d_sel, sizeof(int)*selected.size()));
        V11_CUDA_CHECK(cudaMemcpy(d_sel, selected.data(), sizeof(int)*selected.size(), cudaMemcpyHostToDevice));
        V11_CUDA_CHECK(cudaMalloc(&dVsel_out, sizeof(cuDoubleComplex) * static_cast<std::size_t>(N) * static_cast<std::size_t>(vdim_out)));
        const int vblocks = (N*vdim_out + threads - 1) / threads;
        gather_v32y_evecs_to_complex_kernel<<<vblocks, threads>>>(dVsel_out, d_P, d_sel, N, vdim_out, chop_tol);
        V11_CUDA_CHECK(cudaGetLastError());
        V11_CUDA_CHECK(cudaDeviceSynchronize());

        if (debug == 'y') {
            std::cout << "[projection exact-v32y GPU] N=" << N
                      << " vdim=" << vdim_out
                      << " irrep=" << irrep
                      << " LG=" << nrot
                      << " eval_min=" << evals.front()
                      << " eval_max=" << evals.back() << "\n";
        }

        V11_CUSOLVER_CHECK(cusolverDnDestroy(solver)); solver = nullptr;
    } catch (...) {
        if (solver) cusolverDnDestroy(solver);
        cudaFree(d_nx); cudaFree(d_ny); cudaFree(d_nz); cudaFree(d_ell); cudaFree(d_m); cudaFree(d_flavor);
        cudaFree(d_R); cudaFree(d_sel); cudaFree(d_info); cudaFree(d_chi); cudaFree(d_P); cudaFree(d_W); cudaFree(d_work);
        if (dVsel_out) { cudaFree(dVsel_out); dVsel_out=nullptr; }
        vdim_out = 0;
        throw;
    }

    if (solver) cusolverDnDestroy(solver);
    cudaFree(d_nx); cudaFree(d_ny); cudaFree(d_nz); cudaFree(d_ell); cudaFree(d_m); cudaFree(d_flavor);
    cudaFree(d_R); cudaFree(d_sel); cudaFree(d_info); cudaFree(d_chi); cudaFree(d_P); cudaFree(d_W); cudaFree(d_work);
}

} // namespace v32zt_vsel_gpu

#endif // PROJECTIONS_GPU_V32Y_EXACT_DEVICE_CUH
