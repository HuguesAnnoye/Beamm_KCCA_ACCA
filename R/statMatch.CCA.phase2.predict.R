#' @title statMatch.CCA.phase2.predict
#' @description Prediction of phase 2 of the CCA matching process
#' @param CV_X_A a matrix with the canonical variables for data set A (donnor)
#' @param CV_X_B a matrix with the canonical variables for data set B (receiver)
#' @param h a positive double, value of the bandwidth parameters.
#' @param d numeric value, number of canonical variable to compute
#' @param mtx.rec.CV.raw numeric matrix, common variables in the receiver data set (categorical variables transformed into dummies).
#' @param mtx.don.CV.raw numeric matrix, common variables in the donor data set (categorical variables transformed into dummies).
#' @param mtx.don.NCV.raw numeric matrix, non-common variables in the donor data set (categorical variables transformed into dummies).
#' @param don.weights numeric vector, sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.CCA.options()} which contains all the options.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param keep_unconstrained a logical, whether or not the predictions are adjusted with the zero constraints.
#' @return RETURN
#' @details DETAILS
#' @importFrom magrittr '%>%'
#' @importFrom stats cancor
#' @noRd

statMatch.CCA.phase2.predict <- function(CV_X_A, CV_X_B, h, d, mtx.rec.CV.raw, mtx.don.CV.raw,
                                          mtx.don.NCV.raw, don.weights, opts, zero.constraints, keep_unconstrained = F) {

  # Go back to the original scale of the variables
  if (opts$scaling == "z-score") {
    mtx.rec.CV.raw <- scale_mtx_inv(mtx.rec.CV.raw,
                                    center = attr(mtx.rec.CV.raw, "scaled:center"),
                                    scale = attr(mtx.rec.CV.raw, "scaled:scale"))
    mtx.don.CV.raw <- scale_mtx_inv(mtx.don.CV.raw,
                                    center = attr(mtx.don.CV.raw, "scaled:center"),
                                    scale = attr(mtx.don.CV.raw, "scaled:scale"))
    mtx.don.NCV.raw <- scale_mtx_inv(mtx.don.NCV.raw,
                                     center = attr(mtx.don.NCV.raw, "scaled:center"),
                                     scale = attr(mtx.don.NCV.raw, "scaled:scale"))
  } else if (opts$scaling == "min-max"  & length(zero.constraints) > 0) {
    mtx.rec.CV.raw <- normalize_mtx_inv(mtx.rec.CV.raw,
                                        min = attr(mtx.rec.CV.raw, "min"),
                                        max = attr(mtx.rec.CV.raw, "max"))
    mtx.don.CV.raw <- normalize_mtx_inv(mtx.don.CV.raw,
                                        min = attr(mtx.don.CV.raw, "min"),
                                        max = attr(mtx.don.CV.raw, "max"))
    mtx.don.NCV.raw <- normalize_mtx_inv(mtx.don.NCV.raw,
                                         min = attr(mtx.don.NCV.raw, "min"),
                                         max = attr(mtx.don.NCV.raw, "max"))
  }
  # Compute predictions from CCA
  names.NCV.don <- colnames(mtx.don.NCV.raw)
  is.constrained <- names.NCV.don %in% zero.constraints
  PRED <- matrix(NA_real_, nrow(mtx.rec.CV.raw), ncol(mtx.don.NCV.raw))
  for (i in seq_len(nrow(mtx.rec.CV.raw))) {
    # Compute the euclidean distance to the individuals in the donor data set
    distances <- compute_eucl_dist_vec_mat_rcpp(CV_X_B[i,], CV_X_A)
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
        if (mean(prob_i * !constraint)==0) PRED[i,j] <- 0
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
