#' @title cramer von minses discrete
#' @description cramer von minses between the multivariate and the indepedence case
#' @param DATA a data frame with observations.
#' @param wt a numeric vector, sample weights for ACTUAL.
#' @param names a character vector with the names of the non-common variables used for the prediction.
#' @param type a string, type of error measure, (either \code{mean}, \code{sum} or \code{sumst}).
#' @param parallel a logical, whether or not use parallel computing.
#' @param nnodes an integer, number of parallel session (CPUs).
#' @return A list,with the values of the objective function.
#' @importFrom bivariate ebvcdf
#' @importFrom foreach %dopar%
#' @export

compute_cramer_von_mises_for_independence <- function(DATA, wt = NULL, names=NULL, type = c("mean", "sum", "sumst"), parallel = FALSE, nnodes = 10) {
  if (is.null(names)) names <- colnames(DATA)
  if (!is.null(wt) & (length(wt) != nrow(DATA))) stop("Problem in weights for DATAL.")
  ACTUAL <- df2mtx(DATA[, names])
  names.ACTUAL <- colnames(ACTUAL)
  n_var <- length(names.ACTUAL)
  #matA <- NA

  if (is.null(wt)) {
    wt <- rep(1,nrow(ACTUAL))
  }
  Fn <- ebvcdfwall(ACTUAL, wt = wt)
  Gn <- ebvcdfwind(ACTUAL, wt = wt)
  if (type == "mean") {
     wt2 <- wt / sum(wt)
  } else if (type == "sum") {
    wt2 <- wt
  } else if (type == "sumst") {
    wt2 <- wt / sum(wt) * nrow(ACTUAL)  }
  matA <- sum(wt2*(Fn(ACTUAL) - Gn(ACTUAL))^2)
  #matB <- sum(wt2_PRED*(Fn(PRED) - Gn(PRED))^2)
  output <- list(matA = matA)
  return(output)
}

