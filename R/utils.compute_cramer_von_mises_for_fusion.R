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

compute_cramer_von_mises_for_fusion <- function(ACTUAL, PRED, wt = NULL, wt_PRED=NULL, names.NCV=NULL, type = c("mean", "sum", "sumst"), parallel = FALSE, nnodes = 10) {
  if (is.null(names.NCV)) names.NCV <- colnames(ACTUAL)
  if (!is.null(wt) & (length(wt) != nrow(ACTUAL))) stop("Problem in weights for ACTUAL.")
  if (!is.null(wt_PRED) & (length(wt_PRED) != nrow(PRED))) stop("Problem in weights for PRED.")
  ACTUAL <- df2mtx(ACTUAL[, names.NCV])
  PRED <- df2mtx(PRED[, names.NCV])
  names.ACTUAL <- colnames(ACTUAL)
  names.PRED <- colnames(PRED)
  if (length(names.ACTUAL) != length(names.PRED)) stop("different number of variables in ACTUAL and PRED.")
  if (any(sort(names.ACTUAL) != sort(names.PRED))) stop("diffent variable names in ACTUAL and PRED.")
  n_var <- length(names.ACTUAL)

  if (!parallel) {
    matA <- matrix(NA, nrow = n_var, ncol = n_var)
    matB <- matrix(NA, nrow = n_var, ncol = n_var)
    colnames(matA) <- names.ACTUAL
    rownames(matA) <- names.ACTUAL
    colnames(matB) <- names.ACTUAL
    rownames(matB) <- names.ACTUAL
    for (var in 1:n_var) {
      for (var2 in var:n_var) {
        if (is.null(wt) & is.null(wt_PRED)) {
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
          if (is.null(wt)) wt <- rep(1,nrow(ACTUAL))
          if (is.null(wt_PRED)) wt_PRED <- rep(1,nrow(PRED))
          Fn <- ebvcdfw(ACTUAL[, var], ACTUAL[, var2], wt = wt)
          Gn <- ebvcdfw(PRED[, var], PRED[, var2], wt = wt_PRED)
          # Fn_2 <- ebvcdf(ACTUAL[, var], ACTUAL[, var2])
          # Gn_2 <- ebvcdf(PRED[, var], PRED[, var2])
          if (type == "mean") {
            wt2 <- wt / sum(wt)
            wt2_PRED <- wt_PRED / sum(wt_PRED)
          } else if (type == "sum") {
            wt2 <- wt
            wt2_PRED <- wt_PRED
          } else if (type == "sumst") {
            wt2 <- wt / sum(wt) * nrow(ACTUAL)
            wt2_PRED <- wt_PRED / sum(wt_PRED) * nrow(PRED)
          }
          tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2])
          tmpPRED <- cbind(PRED[,var], PRED[,var2])
          matA[var, var2] <- sum(wt2*(Fn(tmpACTUAL[,1], tmpACTUAL[,2]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2]))^2)
          matB[var, var2] <- sum(wt2_PRED*(Fn(tmpPRED[,1], tmpPRED[,2]) - Gn(tmpPRED[,1], tmpPRED[,2]))^2)
          matA[var2, var] <- matA[var, var2]
          matB[var2, var] <- matB[var, var2]
        }
      }
    }
    matA_final <- matA
    matB_final <- matB
  } else {
    # Initialize cluster
    nnodes <- min(nnodes, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(nnodes, outfile = "")
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    coord <- expand.grid(
      row = 1:n_var,
      col = 1:n_var
    )

    i <- NULL
    matA <- foreach::foreach(i = seq_len(nrow(coord)), .combine = rbind) %dopar% {
      var <- coord$row[i]
      var2 <- coord$col[i]
      # message("matA ",i, " of ",nrow(coord))
      if (var > var2) matA_ij <- NA
      else {
        if (is.null(wt) & is.null(wt_PRED)) {
          Fn <- ebvcdf(ACTUAL[, var], ACTUAL[, var2])
          Gn <- ebvcdf(PRED[, var], PRED[, var2])
          tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2])
          tmpPRED <- cbind(PRED[,var], PRED[,var2])
          if (type == "mean") {
            matA_ij <- mean((Fn(tmpACTUAL[,1], tmpACTUAL[,2]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2]))^2)
          } else {
            matA_ij <- sum((Fn(tmpACTUAL[,1], tmpACTUAL[,2]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2]))^2)
          }
        } else {
          if (is.null(wt)) wt <- rep(1,nrow(ACTUAL))
          if (is.null(wt_PRED)) wt_PRED <- rep(1,nrow(PRED))
          Fn <- ebvcdfw(ACTUAL[, var], ACTUAL[, var2], wt = wt)
          Gn <- ebvcdfw(PRED[, var], PRED[, var2], wt = wt_PRED)
          if (type == "mean") {
            wt2 <- wt / sum(wt)
            wt2_PRED <- wt_PRED / sum(wt_PRED)
          } else if (type == "sum") {
            wt2 <- wt
            wt2_PRED <- wt_PRED
          } else if (type == "sumst") {
            wt2 <- wt / sum(wt) * nrow(ACTUAL)
            wt2_PRED <- wt_PRED / sum(wt_PRED) * nrow(PRED)
          }
          tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2])
          tmpPRED <- cbind(PRED[,var], PRED[,var2])
          matA_ij <- sum(wt2*(Fn(tmpACTUAL[,1], tmpACTUAL[,2]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2]))^2)
        }
      }
      c(var, var2, matA_ij)
    }
    matB <- foreach::foreach(i = seq_len(nrow(coord)), .combine = rbind) %dopar% {
      var <- coord$row[i]
      var2 <- coord$col[i]
      # message("matB ",i, " of ",nrow(coord))
      if (var > var2) matB_ij <- NA
      else {
        if (is.null(wt) & is.null(wt_PRED)) {
          Fn <- ebvcdf(ACTUAL[, var], ACTUAL[, var2])
          Gn <- ebvcdf(PRED[, var], PRED[, var2])
          tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2])
          tmpPRED <- cbind(PRED[,var], PRED[,var2])
          if (type == "mean") {
            matB_ij <- mean((Fn(tmpPRED[,1], tmpPRED[,2]) - Gn(tmpPRED[,1], tmpPRED[,2]))^2)
          } else {
            matB_ij <- sum((Fn(tmpPRED[,1], tmpPRED[,2]) - Gn(tmpPRED[,1], tmpPRED[,2]))^2)
          }
        } else {
          if (is.null(wt)) wt <- rep(1,nrow(ACTUAL))
          if (is.null(wt_PRED)) wt_PRED <- rep(1,nrow(PRED))
          Fn <- ebvcdfw(ACTUAL[, var], ACTUAL[, var2], wt = wt)
          Gn <- ebvcdfw(PRED[, var], PRED[, var2], wt = wt_PRED)
          if (type == "mean") {
            wt2 <- wt / sum(wt)
            wt2_PRED <- wt_PRED / sum(wt_PRED)
          } else if (type == "sum") {
            wt2 <- wt
            wt2_PRED <- wt_PRED
          } else if (type == "sumst") {
            wt2 <- wt / sum(wt) * nrow(ACTUAL)
            wt2_PRED <- wt_PRED / sum(wt_PRED) * nrow(PRED)
          }
          tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2])
          tmpPRED <- cbind(PRED[,var], PRED[,var2])
          matB_ij <- sum(wt2_PRED*(Fn(tmpPRED[,1], tmpPRED[,2]) - Gn(tmpPRED[,1], tmpPRED[,2]))^2)
        }
      }
      c(var, var2, matB_ij)
    }
    matA_final <- matrix(matA[,3], nrow = max(matA[,1]), ncol = max(matA[,2]))
    matB_final <- matrix(matB[,3], nrow = max(matB[,1]), ncol = max(matB[,2]))
  }
  colnames(matA_final) <- colnames(matB_final) <- names.ACTUAL
  rownames(matA_final) <- rownames(matB_final) <- names.ACTUAL
  output <- list(matA = matA_final,
                 matB = matB_final)
  return(output)
}
