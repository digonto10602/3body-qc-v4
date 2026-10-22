#include <cstdio>
#include <string>
#include "fast_fclass.hpp"
using namespace qc3_fast;
int main(int argc, char** argv) {
  // args: nmax alpha then triples
  int nmax = std::atoi(argv[1]); double alpha = std::atof(argv[2]);
  long long totb = 0, totc = 0;
  for (int i = 3; i + 2 < argc; i += 3) {
    int c[3] = {std::atoi(argv[i]), std::atoi(argv[i+1]), std::atoi(argv[i+2])};
    auto t = build_ftable(c, nmax, alpha);
    long long b = (long long)t->cl.size() * sizeof(FClass) + (long long)t->grp.size() * sizeof(FGroup) + sizeof(FTable);
    std::printf("%d %d %d classes %zu groups %zu points %lld bytes %lld\n", c[0], c[1], c[2], t->cl.size(), t->grp.size(), t->total_points, b);
    totb += b; totc += (long long)t->cl.size();
  }
  std::printf("TOTAL classes %lld bytes %lld sizeof(FClass)=%zu sizeof(FGroup)=%zu\n", totc, totb, sizeof(FClass), sizeof(FGroup));
}
