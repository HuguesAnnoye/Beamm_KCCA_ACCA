#' @title compute_mtx_compatibilities
#' @description Compute matrix of compatible observations between df.don and df.rec
#' @param df.rec a data.frame, the recipient data set.
#' @param df.don a data.frame, the donor data set.
#' @param names.CV.cat a character vector, the names of the categorical common variables on which the compatibility matrix is computed.
#' @param print.info a logical, whether information about the compatibilities are printed to the terminal.
#' @param index a logical, whether only the positions of the compatible individuals in the two data.frames are provided.
#' @return Either a matrix of dimension \code{nrow(df.don)} times \code{nrow(df.rec)} or a matrix with two columns with the positions of the ones in the compatibility matrix (to save memory).
#' @noRd

compute_mtx_compatibilities <- function(df.rec, df.don, names.CV.cat, print.info = FALSE, index = FALSE) {

  mtx.rec.CV.cat <- df2mtx(df.rec %>% select(all_of(names.CV.cat)))
  mtx.don.CV.cat <- df2mtx(df.don %>% select(all_of(names.CV.cat)))
  mtx.n_diff_CV <- compute_mtx_n_diff_CV_rcpp(mtx.don.CV.cat, mtx.rec.CV.cat) / 2 # Divide because there is one dummy for each level !!!
  res <- compute_mtx_compatibilities_rcpp(mtx.n_diff_CV)

  if (print.info) {
    message("\n", res$N.df.rec - res$full.obs.avail, " / ", res$N.df.rec, " observations in df.rec without an exact matching of the categorical common-variables.")
    message("Up to ", res$lowest.N.CV, " variables were discarded to find compatible observations.")
    message("Observations in df.don without a compatible observation in df.rec = ", sum(rowSums(res$comp.mtx) == 0L), " / ", res$N.df.don)
    message("Observations in df.rec without a compatible observation in df.don = ", sum(colSums(res$comp.mtx) == 0L), " / ", res$N.df.rec)
  }

  if (index) res$comp.mtx <- which(res$comp.mtx == TRUE, arr.ind = T)

  return(res$comp.mtx)
}
