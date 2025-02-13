#' @title statMatch.ACCA.phase2.predict
#' @description Prediction of phase 2 of the ACCA matching process
#' @param cca_model cca + mtx
#' @param h a positive double, value of the bandwidth parameters.
#' @param d numeric value, number of canonical variable to compute
#' @param don.weights numeric vector, sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.ACCA.options()} which contains all the options.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param keep_unconstrained a logical, whether or not the predictions are adjusted with the zero constraints.
#' @return RETURN
#' @details DETAILS
#' @importFrom magrittr '%>%'
#' @importFrom keras keras_model
#' @importFrom stats cancor
#' @importFrom pracma Rank
#' @noRd
#'
statMatch.ACCA.phase2.predict <- function(cca_model, h, d, don.weights, opts, zero.constraints, keep_unconstrained = F) {
  can_var_CCA_don <- cca_model$can_var_CCA_don
  can_var_CCA_rec <- cca_model$can_var_CCA_rec
  mtx.rec.CV.raw <- cca_model$mtx.rec.CV.raw
  mtx.don.CV.raw <- cca_model$mtx.don.CV.raw
  mtx.don.NCV.raw <- cca_model$mtx.don.NCV.raw
  # Compute predictions from CCA
  names.NCV.don <- colnames(mtx.don.NCV.raw)
  is.constrained <- names.NCV.don %in% zero.constraints
  PRED <- matrix(NA_real_, nrow(mtx.rec.CV.raw), ncol(mtx.don.NCV.raw))
  for (i in seq_len(nrow(mtx.rec.CV.raw))) {
    # Compute the euclidean distance to the individuals in the donor data set
    distances <- compute_eucl_dist_vec_mat_rcpp(can_var_CCA_rec[i,], can_var_CCA_don)
    # Compute the gaussian kernel
    h_min <- min(distances)
    h_final <- pmax(h, h_min)
    if (h_final == 0) h_final <- 1
    kern_dist <-  (sqrt(2*pi))^(d) * exp(-((distances/h_final)^2)/2)
    prob_i <- don.weights * kern_dist
    for (j in seq_len(ncol(PRED))) {
      if (is.constrained[j]) {
        constraint <- round(mtx.don.CV.raw[, paste0("ZC_", names.NCV.don[j], "...1")])  # <----------- Ensure that the zero constraints are equal 0 or 1 !!!
        PRED[i,j] <- wtd_mean_rcpp(mtx.don.NCV.raw[,j], prob_i * !constraint)
      } else {
        PRED[i,j] <- wtd_mean_rcpp(mtx.don.NCV.raw[,j], prob_i)
      }
    }
  }
  # Apply the zero constraints (if required)
  if (!keep_unconstrained) {
    # Impose the constraint
    for (j in seq_len(length(names.NCV.don))) {
      if (is.constrained[j]) {
        constraint <- round(mtx.rec.CV.raw[, paste0("ZC_", names.NCV.don[j], "...1")])  # <----------- Ensure that the zero constraints are equal 0 or 1 !!!
        PRED[constraint == 1, j] <- 0
      }
    }
  }
  return(PRED)
}
