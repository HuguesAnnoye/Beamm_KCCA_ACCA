#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat compute_omega_kcca_rcpp(const arma::mat& uX,
                                  const arma::mat& uY,
                                  const arma::vec& hmin) {
  int n_row_uX = uX.n_rows;
  int n_row_uY = uY.n_rows;
  int n_col_uX = uX.n_cols;

  arma::mat W(n_row_uX, n_row_uY);
    for (int i = 0; i < n_row_uX; i++) {
      for (int j = 0; j < n_row_uY; j++) {
        // Note that in Rcpp rows and columns start at 0 !
        //W(i,j) = 1 / (2 * datum::pi) *
        //  exp((-std::pow(((uX(i,0) - uY(j,0)) / hmin(j)), 2.0)) / 2.0) *
        //  exp((-std::pow(((uX(i,1) - uY(j,1)) / hmin(j)), 2.0)) / 2.0);
        for (int k = 0; k < n_col_uX; k++) {
          W(i,j) += (-std::pow(((uX(i,k) - uY(j,k)) / hmin(j)), 2.0));
        }
        W(i,j) = 1 / (std::pow(std::sqrt((2 * datum::pi)),n_col_uX)) *
          exp(W(i,j) / 2.0);
      }
    }
    return W;
}
