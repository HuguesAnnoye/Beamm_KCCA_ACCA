#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat compute_hmin_kcca_rcpp(const arma::mat& uX,
                                 const arma::mat& uY,
                                 const double& h) {
  int n_row_uX = uX.n_rows;
  int n_row_uY = uY.n_rows;
  int n_col_uX = uX.n_cols;
  arma::mat hmin(n_row_uX, n_row_uY);
  arma::vec hmin2(n_row_uY);
  for (int i = 0; i < n_row_uX; i++) {
    for (int j = 0; j < n_row_uY; j++) {
      // Note that in Rcpp rows and columns start at 0 !
      //hmin(i,j) = sqrt(std::pow(uX(i,0) - uY(j,0), 2.0) + std::pow(uX(i,1) - uY(j,1), 2.0));
      for (int k = 0; k < n_col_uX; k++) {
        hmin(i,j) += std::pow(uX(i,k) - uY(j,k), 2.0);
      }
      hmin(i,j) = std::sqrt(hmin(i,j));
    }
  }
  for (int j = 0; j < n_row_uY; j++) {
    hmin2(j) = min(hmin.col(j));
    if (hmin2(j) < h ) hmin2(j) = h;
  }
  return hmin2;
}

