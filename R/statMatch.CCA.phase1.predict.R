#' @title statMatch.CCA.phase1.predict
#' @description Prediction of phase 1 of the ACCA matching process
#' @param CV_X_A a matrix with the canonical variables for data set A (donnor)
#' @param CV_X_B a matrix with the canonical variables for data set B (receiver)
#' @param h a positive double, value of the bandwidth parameters.
#' @param d numeric value, number of canonical variable to compute
#' @param Y numeric matrix, non-common variables in the donor data set (categorical variables transformed into dummies).
#' @param comp.mtx integer matrix of size N.don-by-N.rec containg a value 1L if an individual in the donor data set is compatible with another one in the receiver data set.
#' @param don.weights numeric vector, sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.ACCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#' @importFrom magrittr '%>%'
#' @importFrom stats cancor
#' @noRd

statMatch.CCA.phase1.predict <- function(CV_X_A, CV_X_B, h, d, Y, comp.mtx, don.weights, opts) {


  # Based on the predictions, pick up at random among the compatible individuals in the donor
  N.rec <- ncol(comp.mtx)
  N.don <- nrow(comp.mtx)
  idx <- 1:N.don
  mtx_PRED <- matrix(NA_real_, N.rec, ncol(Y))
  for (i in seq_len(N.rec)) {
    prob_i <- don.weights
    idx.don.comp <- idx[comp.mtx[,i] == 1L]
    # Compute the euclidean distance to the compatible individuals in the donor data set
    distances <- compute_eucl_dist_vec_mat_rcpp(CV_X_B[i,], CV_X_A[idx.don.comp, , drop = F])
    # Compute the gaussian kernel
    h_min <- min(distances)
    h_final <- pmax(h, h_min)
    if (h_final == 0) h_final <- 1
    kern_dist <-  (sqrt(2*pi))^(d) * exp(-((distances/h_final)^2)/2)
    # Final weights used to pick up at random
    prob_i[idx.don.comp] <- kern_dist * prob_i[idx.don.comp]
    prob_i[-idx.don.comp] <- 0
    # Pick up an individual at random
    mtx_PRED[i,] <- Y[sample(N.don, 1, prob = prob_i),]
  }

  if (opts$scaling == "z-score") {
    mtx_PRED <- scale_mtx_inv(mtx_PRED,
                              center = attr(Y, "scaled:center"),
                              scale = attr(Y, "scaled:scale"))
  } else if (opts$scaling == "min-max") {
    mtx_PRED <- normalize_mtx_inv(mtx_PRED,
                                  min = attr(Y, "min"),
                                  max = attr(Y, "max"))
  }

  return(mtx_PRED)
}
