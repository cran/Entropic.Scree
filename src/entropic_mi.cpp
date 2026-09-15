#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

// Conditionally include OpenMP if the compiler supports it
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(Rcpp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fast_parallel_MI(IntegerMatrix mat, int num_bins, int cores) {
  int n = mat.nrow();
  int p = mat.ncol();
  NumericMatrix MI(p, p);

  std::vector<std::vector<int>> margins(p, std::vector<int>(num_bins, 0));
  std::vector<double> H(p, 0.0);

  for(int j = 0; j < p; ++j) {
    for(int i = 0; i < n; ++i) {
      int val = mat(i, j) - 1;
      if(val >= 0 && val < num_bins) margins[j][val]++;
    }
    double entropy = 0.0;
    for(int b = 0; b < num_bins; ++b) {
      if(margins[j][b] > 0) {
        double prob = (double)margins[j][b] / n;
        entropy -= prob * log(prob);
      }
    }
    H[j] = entropy;
    MI(j, j) = entropy;
  }

  // Apply the parallel pragma only if OpenMP is available
#ifdef _OPENMP
#pragma omp parallel num_threads(cores)
#endif
{
  std::vector<int> joint(num_bins * num_bins, 0);

  // Apply the dynamic schedule only if OpenMP is available
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for(int i = 0; i < p; ++i) {
    for(int j = i + 1; j < p; ++j) {

      std::fill(joint.begin(), joint.end(), 0);

      for(int row = 0; row < n; ++row) {
        int val1 = mat(row, i) - 1;
        int val2 = mat(row, j) - 1;
        if(val1 >= 0 && val1 < num_bins && val2 >= 0 && val2 < num_bins) {
          joint[val1 * num_bins + val2]++;
        }
      }

      double joint_entropy = 0.0;
      for(int k = 0; k < num_bins * num_bins; ++k) {
        if(joint[k] > 0) {
          double prob = (double)joint[k] / n;
          joint_entropy -= prob * log(prob);
        }
      }

      double mi_val = H[i] + H[j] - joint_entropy;
      if (mi_val < 0) mi_val = 0;

      MI(i, j) = mi_val;
      MI(j, i) = mi_val;
    }
  }
}
return MI;
}
