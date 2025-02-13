#' @title Predict for CCA and KCCA
#' @description Do prediction
#' @param W_keep a matrix containing the kernel product BEAMM.KCCAACCA.KCCA.predict3
#' @param Z a matrix with the non-common variables in data set 1
#' @param d a positive integer, number of latent variable used in CCA.
#' @param weights a vector of individual weights
#' @param kernel_predict Type of kernel use for prediction
#' @param rot a logical, if TRUE the bandwidth is h multiply by the variance
#' @param print.details if TRUE print the details
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param matrix.tot.possibilitities A matrix where compatibilities are encodded
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @noRd

# attention il faut tenir compte du lenght variable =1
statMatch.KCCA.phase2.predict.prediction <- function(W_keep,
                                                     Z,
                                                     d = 1,
                                                     kernel_predict = c("gauss", "unif", "epan", "dist", "alea"),
                                                     rot = FALSE,
                                                     matrix.tot.possibilitities = NULL,
                                                     names.NCV = NULL,
                                                     weights = NULL,
                                                     print.details = FALSE) {
  if (is.null(matrix.tot.possibilitities)) {
    matrix.tot.possibilitities <- matrix(1, nrow = nrow(W_keep), ncol = ncol(W_keep))
  }
  W2 <- matrix.tot.possibilitities * W_keep
  matrix.tot.possibilitities <- NULL
  W_keep <- NULL
  Wcs <- colSums(W2)
  Wcs[Wcs == 0] <- 1
  W1 <- t(W2) / Wcs
  W2 <- NULL
  Wcs <- NULL
  if (kernel_predict == "gauss") {
    if (d >= 1) {
      data_testchap <- W1 %*% Z
      W1 <- NULL
      mtx_PRED <- data_testchap[, names.NCV, drop = F]
    }
  } else if (kernel_predict == "alea") {
    if (d >= 1) {
      data_testchap <- matrix(NA, nrow = ncol(W2), ncol = length(names.NCV))
      colnames(data_testchap) <- names.NCV
      for (j in 1:ncol(W2)) {
        iraw <- sample(1:nrow(W2), 1, prob = W1[j, ])
        data_testchap[j, ] <- Z[iraw, names.NCV, drop = F]
      }
      W1 <- NULL
      mtx_PRED <- data_testchap
    }
  } else if (kernel_predict == "epan") {
    if (d == 1) {
      data_testchap <- W1 %*% Z
      W1 <- NULL
      mtx_PRED <- data_testchap[, names.NCV, drop = F]
    } else if (d == 2) {
      data_testchap <- W1 %*% Z
      W1 <- NULL
      mtx_PRED <- data_testchap[, names.NCV, drop = F]
    }
  }
  return(mtx_PRED)
}
