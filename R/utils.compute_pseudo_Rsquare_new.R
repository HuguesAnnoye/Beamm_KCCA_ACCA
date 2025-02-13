#' @title ccompute_pseudo_Rsquare
#' @description Compute the value of the compute_pseudo_Rsquare ans make the mean or the sum
#' @param ACTUAL a numeric matrix, actual observations.
#' @param PRED a numeric matrix, predictions.
#' @param wt a numeric vector, sample weights.
#' @param type a string, type of error measure, (either \code{mean} or \code{sum}).
#' @return A double, the value of the objective function.
#' @details The standardized version of the objective function, computed using \code{type = "wsRMSE"}, enables us to deal with the fact the continuous variables are measured on different scales.
#'
#' @export
compute_pseudo_Rsquare_new <- function(ACTUAL, PRED, wt, type = c("mean", "sum")) {
  type <- match.arg(type)
  var_ACTUAL <- wtd_col_var_rcpp(ACTUAL, wt)
  length_PRED <- sum(wt) #nrow(PRED)
  error_PRED <- wtd_col_means_rcpp(((PRED - ACTUAL)^2),wt)*(length_PRED/(length_PRED - 1))
  if (any(var_ACTUAL == 0)) {
    var_ACTUAL[var_ACTUAL == 0] <- 1
  }
  if (type == "mean") {
    val.objfun <- mean(1 - (error_PRED / var_ACTUAL))
  }
  else if (type == "sum") {
    val.objfun <- sum(1 - (error_PRED / var_ACTUAL))
  }
  return(val.objfun)
}
