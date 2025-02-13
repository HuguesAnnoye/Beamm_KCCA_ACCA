#' @title statMatch.ACCA.phase1.predict
#' @description Prediction of phase 1 of the ACCA matching process
#' @param cca_model keras model, autoencoder trained on the common variables.
# @param autoencoder_Y keras model, autoencoder trained on the non-common variables.
#' @param h a positive double, value of the bandwidth parameters.
# @param mtx.rec.CV.raw numeric matrix, common variables in the receiver data set (categorical variables transformed into dummies).
# @param mtx.don.CV.raw numeric matrix, common variables in the donor data set (categorical variables transformed into dummies).
#' @param mtx.don.NCV.raw numeric matrix, non-common variables in the donor data set (categorical variables transformed into dummies).
#' @param comp.mtx integer matrix of size N.don-by-N.rec containg a value 1L if an individual in the donor data set is compatible with another one in the receiver data set.
#' @param don.weights numeric vector, sample weights of the donor data set.
#' @param d numeric value, number of canonical variable to compute
#' @param opts a list returned by the function \code{statMatch.ACCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#' @importFrom magrittr '%>%'
#' @importFrom keras keras_model
#' @importFrom stats cancor
#' @importFrom pracma Rank
#' @noRd
#'
statMatch.ACCA.phase1.predict <- function(cca_model, h, d,
                                          mtx.don.NCV.raw, comp.mtx, don.weights, opts) {
  #can_var_CCA_mess <- cca_model$can_var_CCA_mess
  can_var_CCA_don <- cca_model$can_var_CCA_don
  can_var_CCA_rec <- cca_model$can_var_CCA_rec
  # Based on the predictions, pick up at random among the compatible individuals in the donor
  N.rec <- ncol(comp.mtx)
  N.don <- nrow(comp.mtx)
  idx <- 1:N.don
  mtx_PRED <- matrix(NA_real_, N.rec, ncol(mtx.don.NCV.raw))
  for (i in seq_len(N.rec)) {
    prob_i <- don.weights
    idx.don.comp <- idx[comp.mtx[, i] == 1L]
    # Compute the euclidean distance to the compatible individuals in the donor data set
    distances <- compute_eucl_dist_vec_mat_rcpp(can_var_CCA_rec[i, ], can_var_CCA_don[idx.don.comp, , drop = F])
    # Compute the gaussian kernel
    h_min <- min(distances)
    h_final <- pmax(h, h_min)
    if (h_final == 0) h_final <- 1
    kern_dist <- (sqrt(2 * pi))^(d) * exp(-((distances / h_final)^2) / 2)
    # Final weights used to pick up at random
    prob_i[idx.don.comp] <- kern_dist * prob_i[idx.don.comp]
    prob_i[-idx.don.comp] <- 0
    # Pick up an individual at random
    mtx_PRED[i, ] <- mtx.don.NCV.raw[sample(N.don, 1, prob = prob_i), ]
  }

  if (opts$scaling == "z-score") {
    mtx_PRED <- scale_mtx_inv(mtx_PRED,
      center = attr(mtx.don.NCV.raw, "scaled:center"),
      scale = attr(mtx.don.NCV.raw, "scaled:scale")
    )
  } else if (opts$scaling == "min-max") {
    mtx_PRED <- normalize_mtx_inv(mtx_PRED,
      min = attr(mtx.don.NCV.raw, "min"),
      max = attr(mtx.don.NCV.raw, "max")
    )
  }

  return(mtx_PRED)
}
