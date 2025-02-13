#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat compute_omega_kcca_d1_rcpp(const arma::vec& uX,
                                  const arma::vec& uY,
                                  const arma::vec& hmin) {
  int n_row_uX = uX.n_elem;
  int n_row_uY = uY.n_elem;
  arma::mat W(n_row_uX, n_row_uY);
  for (int i = 0; i < n_row_uX; i++) {
    for (int j = 0; j < n_row_uY; j++) {
      // Note that in Rcpp rows and columns start at 0 !
      W(i,j) = 1 / (sqrt(2 * datum::pi)) *
        exp((-std::pow(((uX(i) - uY(j)) / hmin(j)), 2.0)) / 2.0);
    }
  }
  return W;
}
