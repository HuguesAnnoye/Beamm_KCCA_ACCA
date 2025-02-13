#' @title cramer von minsen discrete
#' @description Compute the value of the compute_pseudo_Rsquare ans make the mean or the sum
#' @param ACTUAL a data frame, actual observations.
#' @param PRED a data frame, predictions.
#' @param wt a numeric vector, sample weights for ACTUAL.
#' @param wt_PRED a numeric vector, sample weights for ACTUAL.
#' @param names.NCV a character vector with the names of the non-common variables used for the prediction.
#' @param type a string, type of error measure, (either \code{mean}, \code{sum} or \code{sumst}).
#' @param parallel a logical, whether or not use parallel computing.
#' @param nnodes an integer, number of parallel session (CPUs).
#' @return A list,with the values of the objective function.
#' @importFrom bivariate ebvcdf
#' @importFrom foreach %dopar%
#' @export

compute_cramer_von_mises_for_fusion_all_vec <- function(ACTUAL, PRED, square=TRUE, wt = NULL, wt_PRED=NULL, names.NCV=NULL, parallel = FALSE, nnodes = 10) {
  #type <- match.arg(type)
  if (is.vector(ACTUAL)) {
    if (!is.null(wt) & (length(wt) != length(ACTUAL))) stop("Problem in weights for ACTUAL.")
    if (!is.null(wt_PRED) & (length(wt_PRED) != length(PRED))) stop("Problem in weights for PRED.")
    ACTUAL <- df2mtx(data.frame(ACTUAL))
    PRED <- df2mtx(data.frame(PRED))
  } else {
  if (is.null(names.NCV)) names.NCV <- colnames(ACTUAL)
  if (!is.null(wt) & (length(wt) != nrow(ACTUAL))) stop("Problem in weights for ACTUAL.")
  if (!is.null(wt_PRED) & (length(wt_PRED) != nrow(PRED))) stop("Problem in weights for PRED.")

  ACTUAL <- df2mtx(ACTUAL[, names.NCV])
  PRED <- df2mtx(PRED[, names.NCV])
  }
  names.ACTUAL <- colnames(ACTUAL)
  names.PRED <- colnames(PRED)
  if (length(names.ACTUAL) != length(names.PRED)) stop("different number of variables in ACTUAL and PRED.")
  if (!(is.vector(ACTUAL)|ncol(ACTUAL)==1)) {
    if (any(sort(names.ACTUAL) != sort(names.PRED))) stop("diffent variable names in ACTUAL and PRED.")
  }
  n_var <- length(names.ACTUAL)
  #matA <- NA

  if (is.null(wt) & is.null(wt_PRED)) {
    wt <- rep(1,nrow(ACTUAL))
    wt_PRED <- rep(1,nrow(PRED))
  } else {
   if (is.null(wt)) wt <- rep(1,nrow(ACTUAL))
   if (is.null(wt_PRED)) wt_PRED <- rep(1,nrow(PRED))
  }
  Fn <- ebvcdfwall(ACTUAL, wt = wt)
  Gn <- ebvcdfwall(PRED, wt = wt_PRED)
  #if (type == "mean") {
  #   wt2 <- wt / sum(wt)
  #   wt2_PRED <- wt_PRED / sum(wt_PRED)
  #} else if (type == "sum") {
  #  wt2 <- wt
  #  wt2_PRED <- wt_PRED
  #} else if (type == "sumst") {
  #  wt2 <- wt / sum(wt) * nrow(ACTUAL)
  #  wt2_PRED <- wt_PRED / sum(wt_PRED) * nrow(PRED)
  #}
  if (isTRUE(square)) {
  matA <- (Fn(ACTUAL) - Gn(ACTUAL))^2
  matB <- (Fn(PRED) - Gn(PRED))^2
  output <- list(matA = matA,
                 matB = matB)
  } else {
      matA <- (Fn(ACTUAL) - Gn(ACTUAL))
      matB <- (Fn(PRED) - Gn(PRED))
      output <- list(matA = matA,
                     matB = matB)
  }
  return(output)
}

