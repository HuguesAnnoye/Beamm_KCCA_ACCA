#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat compute_hmin_with_comp_mtx_kcca_d1_rcpp(const arma::vec& uX,
                                               const arma::vec& uY,
                                               const arma::mat& comp_mtx,
                                               const double& h) {
  int n_row_uX = uX.n_elem;
  int n_row_uY = uY.n_elem;
  arma::mat hmin(n_row_uX, n_row_uY);
  arma::vec hmin2(n_row_uY);
  for (int i = 0; i < n_row_uX; i++) {
    for (int j = 0; j < n_row_uY; j++) {
      // Note that in Rcpp rows and columns start at 0 !
      hmin(i,j) = sqrt(std::pow(uX(i) - uY(j), 2.0));
    }
  }
  hmin = comp_mtx % hmin;
  hmin(find(comp_mtx == 0)).fill(datum::inf); // Is it correct ??? not 100% sure !
  for (int j = 0; j < n_row_uY; j++) {
    hmin2(j) = min(hmin.col(j));
    if (hmin2(j) < h ) hmin2(j) = h;
  }
  return hmin2;
}
