#' @title cramer von minsen discrete
#' @description Compute the value of the compute_pseudo_Rsquare ans make the mean or the sum
#' @param ACTUAL a data frame, actual observations.
#' @param PRED a data frame, predictions.
#' @param wt a numeric vector, sample weights.
#' @param names.NCV a character vector with the names of the non-common variables used for the prediction.
#' @param type a string, type of error measure, (either \code{mean}, \code{sum} or \code{sumst}).
#' @return A list,with the values of the objective function.
#' @importFrom bivariate ebvcdf
#' @export

compute_cramer_von_mises <- function(ACTUAL, PRED, wt = NULL, names.NCV=NULL, type = c("mean", "sum", "sumst")) {
  if (is.null(names.NCV)) names.NCV <- colnames(ACTUAL)
  ACTUAL <- df2mtx(ACTUAL[, names.NCV])
  PRED <- df2mtx(PRED[, names.NCV])
  names.ACTUAL <- colnames(ACTUAL)
  names.PRED <- colnames(PRED)
  if (length(names.ACTUAL) != length(names.PRED)) stop("different number of variables in ACTUAL and PRED.")
  if (any(sort(names.ACTUAL) != sort(names.PRED))) stop("diffent variable names in ACTUAL and PRED.")
  n_var <- length(names.ACTUAL)

  matA <- matrix(NA, nrow = n_var, ncol = n_var)
  matB <- matrix(NA, nrow = n_var, ncol = n_var)
  colnames(matA) <- names.ACTUAL
  rownames(matA) <- names.ACTUAL
  colnames(matB) <- names.ACTUAL
  rownames(matB) <- names.ACTUAL
  for (var in 1:n_var) {
    for (var2 in var:n_var) {
      if (is.null(wt)) {
        Fn <- ebvcdf(ACTUAL[, var], ACTUAL[, var2])
        Gn <- ebvcdf(PRED[, var], PRED[, var2])
        tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2])
        tmpPRED <- cbind(PRED[,var], PRED[,var2])
        if (type == "mean") {
          matA[var, var2] <- mean((Fn(tmpACTUAL[,1], tmpACTUAL[,2]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2]))^2)
          matB[var, var2] <- mean((Fn(tmpPRED[,1], tmpPRED[,2]) - Gn(tmpPRED[,1], tmpPRED[,2]))^2)
        } else {
          matA[var, var2] <- sum((Fn(tmpACTUAL[,1], tmpACTUAL[,2]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2]))^2)
          matB[var, var2] <- sum((Fn(tmpPRED[,1], tmpPRED[,2]) - Gn(tmpPRED[,1], tmpPRED[,2]))^2)
        }
        matA[var2, var] <- matA[var, var2]
        matB[var2, var] <- matB[var, var2]
      } else {
        Fn <- ebvcdfw(ACTUAL[, var], ACTUAL[, var2], wt = wt)
        Gn <- ebvcdfw(PRED[, var], PRED[, var2], wt = wt)
        if (type == "mean") {
          wt2 <- wt / sum(wt)
        } else if (type == "sum") {
          wt2 <- wt
        } else if (type == "sumst") {
          wt2 <- wt / sum(wt) * nrow(ACTUAL)
        }
        tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2])
        tmpPRED <- cbind(PRED[,var], PRED[,var2])
        matA[var, var2] <- sum(wt2*(Fn(tmpACTUAL[,1], tmpACTUAL[,2]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2]))^2)
        matB[var, var2] <- sum(wt2*(Fn(tmpPRED[,1], tmpPRED[,2]) - Gn(tmpPRED[,1], tmpPRED[,2]))^2)
        matA[var2, var] <- matA[var, var2]
        matB[var2, var] <- matB[var, var2]
      }
    }
  }
  sumA <- c(
    sum(matA[upper.tri(matA, diag = FALSE)]),
    sum(matA[upper.tri(matA, diag = TRUE)]),
    sum(matA)
  )
  sumB <- c(
    sum(matB[upper.tri(matB, diag = FALSE)]),
    sum(matB[upper.tri(matB, diag = TRUE)]),
    sum(matB)
  )
  meanA <- c(
    mean(matA[upper.tri(matA, diag = FALSE)]),
    mean(matA[upper.tri(matA, diag = TRUE)]),
    mean(matA)
  )
  meanB <- c(
    mean(matB[upper.tri(matB, diag = FALSE)]),
    mean(matB[upper.tri(matB, diag = TRUE)]),
    mean(matB)
  )
  val <- c(sumA[2],sumB[2])
  output <- list(matA = matA, matB = matB, val=val, sumA = sumA, sumB = sumB, meanA = meanA, meanB = meanB)
}
