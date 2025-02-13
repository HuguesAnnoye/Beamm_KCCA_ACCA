#' @title compute_multi_Rsquare
#' @description Compute the multivariate R square
#' @param ACTUAL a data.frame, actual observations.
#' @param PRED a data.frame, predictions.
#' @param weights a numeric vector, sample weights.
#' @param names.NCV a character vector with the names of the variables.
#' @param type a character string to indicate which multivariate R square mut be calculated \code{"chow"} for the one of Chow, \code{"book"} for the one of Cohen and \code{R} for the one of Jones.
#' @return A double, the value of the multivariate R square
#'
#' @importFrom stats cov.wt
#'
#' @export

compute_multi_Rsquare <- function(PRED, ACTUAL, names.NCV, weights = NULL, type = c("chow", "book", "R", "R_bis")) {
  type <- match.arg(type)
  ACTUAL <- df2mtx(ACTUAL[, names.NCV], n_1_levels = TRUE)
  PRED <- df2mtx(PRED[, names.NCV], n_1_levels = TRUE)
  if (is.null(weights)) weights <- rep(1, nrow(ACTUAL))
  names.ACTUAL <- colnames(ACTUAL)
  names.PRED <- colnames(PRED)
  if (length(names.ACTUAL) != length(names.PRED)) stop("different number of variables in ACTUAL and PRED.")
  if (any(sort(names.ACTUAL) != sort(names.PRED))) stop("diffent variable names in ACTUAL and PRED.")
  if (type == "chow") {
    n_var <- length(names.ACTUAL)
    cov_mtx <- cov.wt(cbind(ACTUAL, PRED), wt = weights)$cov
    multR2 <- det(cov_mtx[1:n_var, (n_var + 1):(2 * n_var)])^2 /
      (det(cov_mtx[1:n_var, 1:n_var]) * det(cov_mtx[(n_var + 1):(2 * n_var), (n_var + 1):(2 * n_var)]))
  } else if (type == "book") {
    det_ACTUAL_PRED <- det(cov.wt(cbind(ACTUAL, PRED), wt = weights, cor = T)$cor)
    det_ACTUAL <- det(cov.wt(ACTUAL, wt = weights, cor = T)$cor)
    det_PRED <- det(cov.wt(PRED, wt = weights, cor = T)$cor)
    multR2 <- 1 - det_ACTUAL_PRED / det_ACTUAL / det_PRED
  } else if (type == "R_bis") {
    mean_ACTUAL <- wtd_col_means_rcpp(ACTUAL, wt = weights)
    SST <- wtd_mean_rcpp(rowSums((t(t(ACTUAL) - mean_ACTUAL))^2), wt = weights)
    SSr <- wtd_mean_rcpp(rowSums((t(t(PRED) - mean_ACTUAL))^2), wt = weights)
    multR2 <- (SSr / SST)
  } else {
    mean_ACTUAL <- wtd_col_means_rcpp(ACTUAL, wt = weights)
    SST <- wtd_mean_rcpp(rowSums((t(t(ACTUAL) - mean_ACTUAL))^2), wt = weights)
    SSE <- wtd_mean_rcpp(rowSums((ACTUAL - PRED)^2), wt = weights)
    multR2 <- 1 - (SSE / SST)
  }
  return(multR2)
}
