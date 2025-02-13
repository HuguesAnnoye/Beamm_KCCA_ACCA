#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector wtd_col_sd_rcpp(const NumericMatrix& mtx, const NumericVector& wt) {
  int n_row = mtx.nrow();
  int n_wt = wt.length();
  if (n_row != n_wt) {
    stop("Number of row of 'mtx' must be equal to length of 'wt'.");
  }
  double wt_sum = sum(wt);
  int n_var = mtx.ncol();
  NumericVector out(n_var);
  for (int i = 0; i < n_var; i++) {
    //out(i) = std::sqrt(sum(pow((mtx(_,i)-(sum(mtx(_,i)*wt) / wt_sum)),2.0)*wt) / (wt_sum-1)) ;
    out(i) = std::sqrt(sum(pow((mtx(_,i)-(sum(mtx(_,i)*wt) / wt_sum)),2.0)*wt) / (wt_sum)) ;
  }
  return out;
}
