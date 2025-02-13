#' @title ccompute_pseudo_Rsquare
#' @description Compute the value of the compute_pseudo_Rsquare ans make the mean or the sum
#' @param ACTUAL a numeric matrix, actual observations.
#' @param PRED a numeric matrix, predictions.
#' @param wt a numeric vector, sample weights.
#' @param type a string, type of error measure, (either \code{mean} or \code{sum}).
#' @param ACTUAL_mean if TRUE use  the actual mean in the numerator if FALSE use the predicted mean
#' @return A double, the value of the objective function.
#' @details The standardized version of the objective function, computed using \code{type = "wsRMSE"}, enables us to deal with the fact the continuous variables are measured on different scales.
#'
#' @export
compute_pseudo_Rsquare <- function(ACTUAL, PRED, wt, type = c("mean", "sum"), ACTUAL_mean = FALSE) {
  type <- match.arg(type)
  if (ACTUAL_mean == TRUE) {
    var_ACTUAL <- wtd_col_var_rcpp(ACTUAL, wt)
    mean_ACTUAL <- wtd_col_means_rcpp(ACTUAL, wt)
    length_PRED <- nrow(PRED)
    #var_PRED <- colSums((t(t(PRED) - mean_ACTUAL))^2)/(length_PRED-1)
    var_PRED <- wtd_col_means_rcpp(((t(t(PRED) - mean_ACTUAL))^2),wt)*(length_PRED/(length_PRED-1))
    if (any(var_ACTUAL == 0)) {
      var_ACTUAL[var_ACTUAL == 0] <- 1
    }
    if (type == "mean") {
      val.objfun <- mean(var_PRED / var_ACTUAL)
    }
    else if (type == "sum") {
      val.objfun <- sum(var_PRED / var_ACTUAL)
    }
    return(val.objfun)
  } else {
    var_ACTUAL <- wtd_col_var_rcpp(ACTUAL, wt)
    var_PRED <- wtd_col_var_rcpp(PRED, wt)
    if (any(var_ACTUAL == 0)) {
      var_ACTUAL[var_ACTUAL == 0] <- 1
    }
    if (type == "mean") {
      val.objfun <- mean(var_PRED / var_ACTUAL)
    }
    else if (type == "sum") {
      val.objfun <- sum(var_PRED / var_ACTUAL)
    }
    return(val.objfun)
  }
}
