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

compute_cramer_von_mises_for_fusion_all <- function(ACTUAL, PRED, wt = NULL, wt_PRED=NULL, names.NCV=NULL, type = c("mean", "sum", "sumst"), parallel = FALSE, nnodes = 10) {
  type <- match.arg(type)
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
  matA <- sum(wt2*(Fn(ACTUAL) - Gn(ACTUAL))^2)
  matB <- sum(wt2_PRED*(Fn(PRED) - Gn(PRED))^2)
  output <- list(matA = matA,
                 matB = matB)
  return(output)
}

#   matA_final <- matA
#   matB_final <- matB
#tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2], ACTUAL[,var3])
#tmpPRED <- cbind(PRED[,var], PRED[,var2], PRED[,var3])
#
# if (!parallel) {
#   matA <- array(NA, dim=c(n_var,n_var,n_var), dimnames = list(names.ACTUAL,names.ACTUAL,names.ACTUAL))
#   matB <- array(NA, dim=c(n_var,n_var,n_var), dimnames = list(names.ACTUAL,names.ACTUAL,names.ACTUAL))
#   for (var in 1:n_var) {
#     for (var2 in var:n_var) {
#       for (var3 in var2:n_var) {
#       if (is.null(wt) & is.null(wt_PRED)) {
#         wt <- rep(1,nrow(ACTUAL))
#         wt_PRED <- rep(1,nrow(PRED))
#         Fn <- ebvcdfw3(ACTUAL[, var], ACTUAL[, var2], ACTUAL[, var3], wt = wt)
#         Gn <- ebvcdfw3(PRED[, var], PRED[, var2], PRED[, var2], wt = wt_PRED)
#         tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2], ACTUAL[,var3])
#         tmpPRED <- cbind(PRED[,var], PRED[,var2], PRED[,var3])
#         if (type == "mean") {
#           matA[var, var2, var3] <- mean((Fn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]))^2)
#           matB[var, var2, var3] <- mean((Fn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]) - Gn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]))^2)
#         } else {
#           matA[var, var2] <- sum((Fn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,2]))^2)
#           matB[var, var2] <- sum((Fn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]) - Gn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]))^2)
#         }
#         matA[var2, var, var3] <- matA[var, var2, var3]
#         matB[var2, var, var3] <- matB[var, var2, var3]
#         matA[var3, var, var2] <- matA[var, var2, var3]
#         matB[var3, var, var2] <- matB[var, var2, var3]
#         matA[var3, var2, var] <- matA[var, var2, var3]
#         matB[var3, var2, var] <- matB[var, var2, var3]
#       } else {
#         if (is.null(wt)) wt <- rep(1,nrow(ACTUAL))
#         if (is.null(wt_PRED)) wt_PRED <- rep(1,nrow(PRED))
#         Fn <- ebvcdfw3(ACTUAL[, var], ACTUAL[, var2], ACTUAL[, var3], wt = wt)
#         Gn <- ebvcdfw3(PRED[, var], PRED[, var2], PRED[, var3], wt = wt_PRED)
#         # Fn_2 <- ebvcdf(ACTUAL[, var], ACTUAL[, var2])
#         # Gn_2 <- ebvcdf(PRED[, var], PRED[, var2])
#         if (type == "mean") {
#           wt2 <- wt / sum(wt)
#           wt2_PRED <- wt_PRED / sum(wt_PRED)
#         } else if (type == "sum") {
#           wt2 <- wt
#           wt2_PRED <- wt_PRED
#         } else if (type == "sumst") {
#           wt2 <- wt / sum(wt) * nrow(ACTUAL)
#           wt2_PRED <- wt_PRED / sum(wt_PRED) * nrow(PRED)
#         }
#         tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2], ACTUAL[,var3])
#         tmpPRED <- cbind(PRED[,var], PRED[,var2], PRED[,var3])
#         matA[var, var2, var3] <- sum(wt2*(Fn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]))^2)
#         matB[var, var2, var3] <- sum(wt2_PRED*(Fn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]) - Gn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]))^2)
#         matA[var2, var, var3] <- matA[var, var2, var3]
#         matB[var2, var, var3] <- matB[var, var2, var3]
#         matA[var3, var, var2] <- matA[var, var2, var3]
#         matB[var3, var, var2] <- matB[var, var2, var3]
#         matA[var3, var2, var] <- matA[var, var2, var3]
#         matB[var3, var2, var] <- matB[var, var2, var3]
#       }
#     }
#     }
#   }
#   matA_final <- matA
#   matB_final <- matB
# } else {
#   # Initialize cluster
#   nnodes <- min(nnodes, parallel::detectCores() - 1)
#   cl <- parallel::makeCluster(nnodes, outfile = "")
#   doParallel::registerDoParallel(cl)
#   on.exit(parallel::stopCluster(cl))
#   coord <- expand.grid(
#     row = 1:n_var,
#     col = 1:n_var,
#     cot = 1:n_var
#   )
#
#   i <- NULL
#   matA <- foreach::foreach(i = seq_len(nrow(coord)), .combine = rbind) %dopar% {
#     var <- coord$row[i]
#     var2 <- coord$col[i]
#     var3 <- coord$cot[i]
#     #message("matA ",i, " of ",nrow(coord))
#     if (var > var2) matA_ij <- NA
#     else if (var2 > var3) matA_ij <- NA
#     else {
#       if (is.null(wt) & is.null(wt_PRED)) {
#         wt <- rep(1,nrow(ACTUAL))
#         wt_PRED <- rep(1,nrow(PRED))
#         Fn <- ebvcdfw3(ACTUAL[, var], ACTUAL[, var2], ACTUAL[, var3], wt = wt)
#         Gn <- ebvcdfw3(PRED[, var], PRED[, var2], PRED[, var3], wt = wt_PRED)
#         tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2], ACTUAL[,var3])
#         #tmpPRED <- cbind(PRED[,var], PRED[,var2], PRED[,var3])
#         if (type == "mean") {
#           matA_ij <- mean((Fn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]))^2)
#         } else {
#           matA_ij <- sum((Fn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]))^2)
#         }
#       } else {
#         if (is.null(wt)) wt <- rep(1,nrow(ACTUAL))
#         if (is.null(wt_PRED)) wt_PRED <- rep(1,nrow(PRED))
#         #message("ebvcdfw3")
#         Fn <- ebvcdfw3(ACTUAL[, var], ACTUAL[, var2], ACTUAL[, var3], wt = wt)
#         Gn <- ebvcdfw3(PRED[, var], PRED[, var2], PRED[, var3], wt = wt_PRED)
#         #message("ebvcdfw3 done")
#         if (type == "mean") {
#           wt2 <- wt / sum(wt)
#           wt2_PRED <- wt_PRED / sum(wt_PRED)
#         } else if (type == "sum") {
#           wt2 <- wt
#           wt2_PRED <- wt_PRED
#         } else if (type == "sumst") {
#           wt2 <- wt / sum(wt) * nrow(ACTUAL)
#           wt2_PRED <- wt_PRED / sum(wt_PRED) * nrow(PRED)
#         }
#         tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2], ACTUAL[,var3])
#        #tmpPRED <- cbind(PRED[,var], PRED[,var2], PRED[,var3])
#         #message(paste0("Mat_A_ij_", i))
#         matA_ij <- sum(wt2*(Fn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]) - Gn(tmpACTUAL[,1], tmpACTUAL[,2], tmpACTUAL[,3]))^2)
#         #message(paste0("Mat_A_ij_done_", i))
#       }
#       #message(paste0("Mat_A_", i))
#     }
#     return(c(var, var2, var3, matA_ij))
#   }
#   message("mat_B")
#   matB <- foreach::foreach(i = seq_len(nrow(coord)), .combine = rbind) %dopar% {
#     var <- coord$row[i]
#     var2 <- coord$col[i]
#     var3 <- coord$cot[i]
#     # message("matB ",i, " of ",nrow(coord))
#     if (var > var2) matB_ij <- NA
#     else if (var2>var3) matB_ij <- NA
#     else {
#       if (is.null(wt) & is.null(wt_PRED)) {
#         wt <- rep(1,nrow(ACTUAL))
#         wt_PRED <- rep(1,nrow(PRED))
#         Fn <- ebvcdfw3(ACTUAL[, var], ACTUAL[, var2], ACTUAL[, var3], wt = wt)
#         Gn <- ebvcdfw3(PRED[, var], PRED[, var2], PRED[, var3], wt = wt_PRED)
#         #tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2], ACTUAL[,var3])
#         tmpPRED <- cbind(PRED[,var], PRED[,var2], PRED[,var3])
#         if (type == "mean") {
#           matB_ij <- mean((Fn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]) - Gn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]))^2)
#         } else {
#           matB_ij <- sum((Fn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]) - Gn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]))^2)
#         }
#       } else {
#         if (is.null(wt)) wt <- rep(1,nrow(ACTUAL))
#         if (is.null(wt_PRED)) wt_PRED <- rep(1,nrow(PRED))
#         Fn <- ebvcdfw3(ACTUAL[, var], ACTUAL[, var2], ACTUAL[, var3], wt = wt)
#         Gn <- ebvcdfw3(PRED[, var], PRED[, var2], PRED[, var3], wt = wt_PRED)
#         if (type == "mean") {
#           wt2 <- wt / sum(wt)
#           wt2_PRED <- wt_PRED / sum(wt_PRED)
#         } else if (type == "sum") {
#           wt2 <- wt
#           wt2_PRED <- wt_PRED
#         } else if (type == "sumst") {
#           wt2 <- wt / sum(wt) * nrow(ACTUAL)
#           wt2_PRED <- wt_PRED / sum(wt_PRED) * nrow(PRED)
#         }
#         #tmpACTUAL <- cbind(ACTUAL[,var], ACTUAL[,var2], ACTUAL[,var3])
#         tmpPRED <- cbind(PRED[,var], PRED[,var2], PRED[,var3])
#         matB_ij <- sum(wt2_PRED*(Fn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]) - Gn(tmpPRED[,1], tmpPRED[,2], tmpPRED[,3]))^2)
#       }
#       #message(paste0("Mat_B_", i))
#     }
#     #message(paste0("Mat_B_", i))
#     return(c(var, var2, var3, matB_ij))
#   }
#   matA_final <- array(matA[,4], dim=c(max(matA[,1]),max(matA[,2]),max(matA[,3])), dimnames = list(names.ACTUAL,names.ACTUAL,names.ACTUAL))
#   matB_final <- array(matB[,4], dim=c(max(matB[,1]),max(matB[,2]),max(matB[,3])), dimnames = list(names.ACTUAL,names.ACTUAL,names.ACTUAL))
#   #matA_final <- matrix(matA[,3], nrow = max(matA[,1]), ncol = max(matA[,2]))
#   #matB_final <- matrix(matB[,3], nrow = max(matB[,1]), ncol = max(matB[,2]))
# }
# #colnames(matA_final) <- colnames(matB_final) <- names.ACTUAL
# #rownames(matA_final) <- rownames(matB_final) <- names.ACTUAL
