#include <cstdio>
#include <string>
#include <random>
#include <set>
#include <algorithm>
#include "fast_fclass.hpp"
using namespace qc3_fast;
typedef long double LD;
struct Mom { LD T0, T1[3], T2[3][3]; long long terms; LD A0,A1,A2; };
// brute force in long double over the whole cube, r^2 <= R2, direct formula (oracle form, nP actual)
Mom brute(const int nP[3], int nmax, double alpha, double x2, double gamma, double xii, double R2) {
  Mom m{}; 
  LD nPsq = nP[0]*nP[0]+nP[1]*nP[1]+nP[2]*nP[2];
  for (int i=-nmax;i<=nmax;++i) for (int j=-nmax;j<=nmax;++j) for (int k=-nmax;k<=nmax;++k) {
    LD r[3] = {(LD)i,(LD)j,(LD)k};
    if (nPsq>0) { LD dot = i*nP[0]+j*nP[1]+k*nP[2]; LD prod1 = dot/nPsq*(1.0L/gamma-1.0L) - (LD)xii/gamma; for(int a=0;a<3;++a) r[a] += nP[a]*prod1; }
    LD r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    if (r2 > (LD)R2) continue;
    LD dd = (LD)x2 - r2; LD w = expl((LD)alpha*dd)/dd;
    m.A0+=fabsl(w); m.A1+=fabsl(w)*sqrtl(r2); m.A2+=fabsl(w)*r2; m.T0 += w; for(int a=0;a<3;++a){ m.T1[a]+=w*r[a]; for(int b=0;b<3;++b) m.T2[a][b]+=w*r[a]*r[b]; }
    m.terms++;
  }
  return m;
}

Mom direct_double(const int nP[3], int nmax, double alpha, double x2, double gamma, double xii, double R2) {
  Mom m{}; double nPsq = nP[0]*nP[0]+nP[1]*nP[1]+nP[2]*nP[2]; bool boosted = nPsq!=0.0;
  double invg_m1 = 1.0/gamma - 1.0, xibygamma = xii/gamma;
  for (int ix=-nmax;ix<=nmax;++ix) for (int iy=-nmax;iy<=nmax;++iy) for (int iz=-nmax;iz<=nmax;++iz) {
    double rx=ix, ry=iy, rz=iz;
    if (boosted) { double dot = ix*(double)nP[0]+iy*(double)nP[1]+iz*(double)nP[2]; double prod1=(dot/nPsq)*invg_m1 - xibygamma; rx=ix+nP[0]*prod1; ry=iy+nP[1]*prod1; rz=iz+nP[2]*prod1; }
    double rsq = rx*rx+ry*ry+rz*rz; if (rsq > R2) continue;
    double dd = x2-rsq; double w = std::exp(alpha*dd)/dd; double r[3]={rx,ry,rz};
    m.T0 += w; for(int a=0;a<3;++a){ m.T1[a]+=w*r[a]; for(int b=0;b<3;++b) m.T2[a][b]+=w*r[a]*r[b]; } m.terms++;
  }
  return m;
}
int main() {
  int nmax = 20; double alpha = 0.5;
  std::mt19937_64 g(20261022);
  struct Case { double x2, gamma, xii; };
  std::vector<Case> cases = {{-2.0,1.3,0.55},{-6.5,2.0,0.5},{0.37,1.1,0.52},{1.9,8.0,0.6},{-0.8,9.98,0.45},{3.2,1.0,0.5},{-40.0,1.2,0.5},{5.5,5.0,0.9}};
  int cs[][3] = {{0,0,0},{1,0,0},{1,1,0},{1,1,1},{2,0,0},{2,1,0},{2,1,1},{3,2,1},{2,2,2},{3,0,0}};
  double worst = 0, worstT = 0, worstD = 0; long long ncmp=0; bool allterms=true;
  for (auto& cs_ : cs) {
    std::shared_ptr<const FTable> tab = ftable_cached(cs_, nmax, alpha);
    for (auto& cc : cases) {
      double R2 = std::max(cc.x2,0.0) + std::log(1e16)/alpha;  // same form as truncation_R2 start; refine
      for (int it=0; it<60; ++it) { double d=R2-cc.x2; double nx = cc.x2 + (std::log(3*R2/d)+std::log(1e16))/alpha; if (std::fabs(nx-R2)<1e-12*R2){R2=nx;break;} R2=nx; }
      R2 *= (1+1e-12);
      FMoments fm = class_moments(*tab, cc.x2, cc.gamma, cc.xii, R2);
      // all distinct signed permutations of c
      std::set<std::array<int,3>> orb;
      int perm[3]={0,1,2};
      do { for (int s=0;s<8;++s) { std::array<int,3> v; for(int a=0;a<3;++a) v[a] = ((s>>a)&1 ? -1:1)*cs_[perm[a]]; orb.insert(v);} } while(std::next_permutation(perm,perm+3));
      for (auto& v : orb) {
        int nPi[3]={v[0],v[1],v[2]};
        int idx[3]={0,1,2}; std::stable_sort(idx,idx+3,[&](int a,int b){return std::abs(nPi[a])>std::abs(nPi[b]);});
        int cv[3], sg[3]; for(int j=0;j<3;++j){cv[j]=std::abs(nPi[idx[j]]); sg[j]=nPi[idx[j]]<0?-1:1;}
        if (cv[0]!=cs_[0]||cv[1]!=cs_[1]||cv[2]!=cs_[2]) { printf("canon mismatch\n"); return 1; }
        double T1[3], T2[3][3];
        for (int j=0;j<3;++j) T1[idx[j]] = sg[j]*fm.T1[j];
        for (int j=0;j<3;++j) for(int k=0;k<3;++k) T2[idx[j]][idx[k]] = sg[j]*sg[k]*fm.T2[j][k];
        Mom bm = brute(nPi, nmax, alpha, cc.x2, cc.gamma, cc.xii, R2);
        if (bm.terms != fm.terms) { allterms=false; printf("TERMS differ c=%d%d%d nP=%d,%d,%d x2=%g: brute %lld class %lld\n",cs_[0],cs_[1],cs_[2],v[0],v[1],v[2],cc.x2,bm.terms,fm.terms); }
        // norm-wise relative errors
        LD n1=0,n2=0, e1=0,e2=0; 
        for(int a=0;a<3;++a){ n1+=bm.T1[a]*bm.T1[a]; e1+=(T1[a]-bm.T1[a])*(T1[a]-bm.T1[a]); for(int b=0;b<3;++b){ n2+=bm.T2[a][b]*bm.T2[a][b]; e2+=(T2[a][b]-bm.T2[a][b])*(T2[a][b]-bm.T2[a][b]); } }
        Mom dm = direct_double(nPi, nmax, alpha, cc.x2, cc.gamma, cc.xii, R2); { LD e0=fabsl(dm.T0-bm.T0)/bm.A0; LD d1=0,d2=0; for(int a=0;a<3;++a){ d1+=(dm.T1[a]-bm.T1[a])*(dm.T1[a]-bm.T1[a]); for(int b=0;b<3;++b) d2+=(dm.T2[a][b]-bm.T2[a][b])*(dm.T2[a][b]-bm.T2[a][b]); } double md=(double)std::max(e0,std::max(sqrtl(d1)/bm.A1,sqrtl(d2)/bm.A2)); worstD=std::max(worstD,md); }
        double r0 = (double)(fabsl(fm.T0-bm.T0)/bm.A0), r1=(double)(sqrtl(e1)/bm.A1), r2=(double)(sqrtl(e2)/bm.A2);
        // scale-aware: error normalised to sum of |terms| unavailable; report norm-relative
        {double mx=std::max(r0,std::max(r1,r2)); if(mx>1e-12) printf("BAD c=%d%d%d nP=%d,%d,%d x2=%g g=%g: r0=%.2e r1=%.2e r2=%.2e T0=%.6Lg\n",cs_[0],cs_[1],cs_[2],v[0],v[1],v[2],cc.x2,cc.gamma,r0,r1,r2,bm.T0); worst = std::max(worst, mx);} ++ncmp;
      }
    }
  }
  printf("direct-double worst = %.3e\n", worstD); printf("compared %lld (nP, case) pairs; worst error/(sum|term|-scale) of (T0,T1,T2) rotated-class vs long-double brute force = %.3e; all term counts equal: %d\n", ncmp, worst, (int)allterms);
}
