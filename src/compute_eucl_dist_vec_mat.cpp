#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_eucl_dist_vec_mat_rcpp(const NumericVector& vec,
                                             const NumericMatrix& mtx) {
  int n_row_mtx = mtx.nrow();
  int n_col_mtx = mtx.ncol();
  NumericVector distances(n_row_mtx);
  distances.fill(0);
  for (int i = 0; i < n_row_mtx; i++) {
    // distances[i] = sqrt(sum(pow(mtx.row(i) - vec, 2.0)));
    for (int j = 0; j < n_col_mtx; j++) {
      distances[i] += std::pow(mtx(i,j) - vec[j], 2.0);
    }
    distances[i] = std::sqrt(distances[i]);
  }
  return distances;
}
