#' @title compute_obj_function
#' @description Compute the value of the objective function (cross-validation), either weighted standardized Root Mean Square Error (\code{wsRMSE}) or weighted Root Mean Square Error (\code{wRMSE}).
#' @param ACTUAL a numeric matrix, actual observations.
#' @param PRED a numeric matrix, predictions.
#' @param wt a numeric vector, sample weights.
#' @param eps a double, machine epsilon.
#' @param type a string, type of error measure, (either \code{wsRMSE} or \code{wRMSE}).
#' @return A double, the value of the objective function.
#' @details The standardized version of the objective function, computed using \code{type = "wsRMSE"}, enables us to deal with the fact the continuous variables are measured on different scales.
#'
#' @export

compute_obj_function <- function(ACTUAL, PRED, wt, type = c("wsRMSE", "wRMSE", "wCE"), eps = .Machine$double.eps) {
  type <- match.arg(type)
  if (type == "wsRMSE") {
    mean_sq_diff <- wtd_col_means_rcpp((ACTUAL - PRED)^2, wt)
    var_ACTUAL <- wtd_col_var_rcpp(ACTUAL, wt)
    if (any(var_ACTUAL == 0))
      var_ACTUAL[var_ACTUAL == 0] <- 1
    val.objfun <- mean(sqrt(mean_sq_diff / var_ACTUAL))
  } else if (type == "wRMSE") {
    mean_sq_diff <- wtd_col_means_rcpp((ACTUAL - PRED)^2, wt)
    val.objfun <- mean(sqrt(mean_sq_diff))
  } else {
    PRED <- pmax(PRED, eps)
    val.objfun <- -wtd_mean_rcpp(rowSums(ACTUAL * log(PRED)), wt)
  }
  return(val.objfun)
}
