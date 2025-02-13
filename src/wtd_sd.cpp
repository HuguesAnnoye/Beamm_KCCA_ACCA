#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double wtd_sd_rcpp(const NumericVector& vec, const NumericVector& wt) {
  int n_vec = vec.length();
  int n_wt = wt.length();
  if (n_vec != n_wt) {
    stop("Length of 'vec' must be equal to length of 'wt'.");
  }
  //double out = std::sqrt(sum(pow((vec-(sum(vec*wt) / sum(wt))),2.0)*wt) / (sum(wt)-1));
  double out = std::sqrt(sum(pow((vec-(sum(vec*wt) / sum(wt))),2.0)*wt) / (sum(wt)));
  return out;
}
