#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List compute_mtx_compatibilities_reg_rcpp(const IntegerMatrix& mtx_n_diff_CV) {
  int n_row_don = mtx_n_diff_CV.nrow();
  int n_row_rec = mtx_n_diff_CV.ncol();
  IntegerMatrix comp_mtx(n_row_don, n_row_rec);
  int full_obs_avail = 0;
  int zero_obs_avail = 0;
  LogicalVector tmp;

  for (int i = 0; i < n_row_rec; i++) {
    tmp = (mtx_n_diff_CV(_,i) == 0);
    if (any(tmp).is_true()) {
      comp_mtx(_,i) = tmp;
      full_obs_avail += 1;
    } else {
      zero_obs_avail += 1;
    }
  }
  return List::create(_["N.df.rec"] = n_row_rec,
                      _["N.df.don"] = n_row_don,
                      _["full.obs.avail"] = full_obs_avail,
                      _["zero.obs.avail"] = zero_obs_avail,
                      _["comp.mtx"] = comp_mtx);
}
