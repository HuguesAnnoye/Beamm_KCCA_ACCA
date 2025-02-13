#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_mtx_n_diff_CV_rcpp(const NumericMatrix& mtx_don,
                                         const NumericMatrix& mtx_rec) {
  int n_row_don = mtx_don.nrow();
  int n_row_rec = mtx_rec.nrow();
  NumericMatrix mtx_n_equal_CV(n_row_don, n_row_rec);
  for (int i = 0; i < n_row_don; i++) {
    for (int j = 0; j < n_row_rec; j++) {
      mtx_n_equal_CV(i,j) = sum(mtx_don.row(i) != mtx_rec.row(j));
    }
  }
  return mtx_n_equal_CV;
}
