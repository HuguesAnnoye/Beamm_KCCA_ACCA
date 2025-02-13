#' @title Fitting for PCA
#' @description Fit with PCA - CCA
#' @param X a matrix with the common variables in data set 1 with n-1 the dummies.
#' @param Y a matrix with the non-common variables in data set 1  with n-1 the dummies.
#' @param X2 a matrix with the common variables in data set 2.
#' @param weights a vector of individual weights.
#' @param d a positive integer, number of latent variable used in CCA.
#' @param lat_X size of the latent space of the common variables.
#' @param lat_Y size of the latent space of the common variables.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom methods is
#' @importFrom stats sd cancor prcomp
#' @noRd

statMatch.PCA.fit <- function(X,
                               Y,
                               X2 = NULL,
                               weights = NULL,
                               d = 1,
                               lat_X = 1,
                               lat_Y = 1) {

  # Weights
  if (is.null(weights)) weights <- rep(1,length = nrow(X))
  Wps <- weights
  sWps <- sum(Wps)
  cca_weights <- sqrt(Wps / sWps)
  # PCA for X
  PCA_X <- prcomp(cca_weights*X,retx = TRUE, center = FALSE, scale. = FALSE,
                  tol = NULL, rank. = lat_X)

  # Prediction
  X_hat <- X %*% PCA_X$rotation
  if (!is.null(X2)) {
    Xt_hat <- X2 %*% PCA_X$rotation
  }

  # PCA for Y
  PCA_Y <- prcomp(cca_weights*Y,retx = TRUE, center = FALSE, scale. = FALSE,
                  tol = NULL, rank. = lat_Y)
  #Y_hat <- PCA_Y$x
  Y_hat <- Y %*% PCA_Y$rotation

  # Centring before CCA
  Xps <- scale_mtx(X_hat, wt = weights)
  Yps <- scale_mtx(Y_hat, wt = weights)
  if (!is.null(X2)) {
  mdX <- attr(Xps, "scaled:center")
  sdX <- attr(Xps, "scaled:scale")
    Xpt <-  t((t(Xt_hat) - mdX) / sdX)
  }
  ## Calculate CCA
  ccas <- cancor(cca_weights * Xps, cca_weights * Yps, xcenter = FALSE, ycenter = FALSE)
  nomx <- names(ccas$xcoef[, 1])
  nlx <- nrow(Xps[, colnames(Xps) %in% nomx])
  ncx <- ncol(Xps[, colnames(Xps) %in% nomx])
  Xps2 <- Xps[, colnames(Xps) %in% nomx]
  Yps <- NULL

  if (d >= 1) {
    # Calculate Canonical variables
    #echap <- cbind((Xps2) %*% ccas$xcoef[, 1], (Xps2) %*% ccas$xcoef[, 2])
    echap <- (Xps2) %*% ccas$xcoef[, 1:d]
    Xps2 <- NULL
    if (!is.null(X2)) {
      Xpt2 <- Xpt[, colnames(Xps) %in% nomx]
      Xpt <- NULL
      Xps <- NULL
      mdx <- NULL
      sdX <- NULL
      echap2 <- (Xpt2) %*% ccas$xcoef[, 1:d]
      Xpt2 <- NULL
      echaplist <- list(CV_X_A = echap, CV_X_B = echap2, CCA = ccas,  PCA_X = PCA_X)
      return(echaplist)
    } else {
      echaplist <- list(CV_X_A = echap, CVec_X_b = ccas$xcoef, CCA = ccas,  PCA_X = PCA_X)
      return(echaplist)
    }
  } else {
    print("d must be a positive integer")
  }
}
