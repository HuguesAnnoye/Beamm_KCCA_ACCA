#' @title Fit for CCA
#' @description Calculate the canonical variables to preform the matching.
#' @param X a matrix with the common variables in data set 1 with n-1 the dummies.
#' @param Y a matrix with the non-common variables in data set 1  with n-1 the dummies.
#' @param X2 a matrix with the common variables in data set 2.
#' @param weights a vector of individual weights.
#' @param d a positive integer, number of latent variable used in CCA.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom stats sd cancor
#'
#' @noRd
statMatch.CCA.fit <- function(X,
                              Y,
                              X2 = NULL,
                              d = 1,
                              weights = NULL) {

  if (is.null(weights)) {

    # Standardized data
    Xps <- scale_mtx(X)
    Yps <- scale_mtx(Y)
    mdX <- attr(Xps, "scaled:center")
    sdX <- attr(Xps, "scaled:scale")
    Y <- NULL

    ## CCA
    ccas <- cancor(Xps, Yps)
    nomy <- names(ccas$ycoef[, 1])
    nomx <- names(ccas$xcoef[, 1])
    nlx <- nrow(Xps[, colnames(Xps) %in% nomx])
    ncx <- ncol(Xps[, colnames(Xps) %in% nomx])
    mxcenter <- matrix(data = ccas$xcenter[names(ccas$xcenter) %in% nomx], nrow = nlx, ncol = ncx, byrow = TRUE)

    Xps2 <- Xps[, colnames(Xps) %in% nomx]
    Yps <- NULL

    if (d == 1) {
      # Calculate Canonical variables
      echap <- (Xps2 - mxcenter) %*% ccas$xcoef[, 1]
      if (!is.null(X2)) {
        Xpt <- t(t(X2) / sdX)
        nlx2 <- nrow(Xpt[, colnames(Xps) %in% nomx])
        mxcenter2 <- matrix(data = ccas$xcenter[names(ccas$xcenter) %in% nomx], nrow = nlx2, ncol = ncx, byrow = TRUE)
        Xpt2 <- Xpt[, colnames(Xps) %in% nomx]
        echap2 <- (Xpt2 - mxcenter2) %*% ccas$xcoef[, 1]
        Xpt2 <- NULL
        mxcenter2 <- NULL
        echaplist <- list(CV_X_A = echap, CV_X_B = echap2, CCA = ccas)
        ccas <- NULL
        echap <- NULL
        echap2 <- NULL
        return(echaplist)
      }
      else {
        echaplist <- list(CV_X_A = echap, Y_A_cv=cbind(Y,echap=echap), CVec_X_B = ccas$xcoef, CCA = ccas)
        echap <- NULL
        ccas <- NULL
        return(echaplist)
      }
    } else if (d >= 2) {
      # Calculate Canonical variables
      echap <- (Xps2 - mxcenter) %*% ccas$xcoef[, 1:d]
      if (!is.null(X2)) {
        Xpt <- t(t(X2) / sdX)
        nlx2 <- nrow(Xpt[, colnames(Xps) %in% nomx])
        mxcenter2 <- matrix(data = ccas$xcenter[names(ccas$xcenter) %in% nomx], nrow = nlx2, ncol = ncx, byrow = TRUE)
        Xpt2 <- Xpt[, colnames(Xps) %in% nomx]
        Xpt <- NULL
        Xps <- NULL
        sdX <- NULL
        echap2 <- (Xpt2 - mxcenter2) %*% ccas$xcoef[, 1:d]
        Xpt2 <- NULL
        mxcenter2 <- NULL
        echaplist <- list(CV_X_A = echap, CV_X_B = echap2, CCA = ccas)
        return(echaplist)
      }
      else {
        echaplist <- list(CV_X_A = echap, Y_A_cv=cbind(Y,echap=echap), CVec_X_B = ccas$xcoef, CCA = ccas)
        return(echaplist)
      }
    }

    else {
      print("d must be a positive integer")
    }
  } else {

    # Weights
    Wps <- weights
    sWps <- sum(Wps)
    cWps <- sqrt(Wps / sWps)
    # cWps <- sqrt(Wps)

    # Standardized data
    Xps <- scale_mtx(X, wt = Wps)
    Yps <- scale_mtx(Y, wt = Wps)
    mdX <- attr(Xps, "scaled:center")
    sdX <- attr(Xps, "scaled:scale")

    ## Calculate CCA
    ccas <- cancor(cWps * Xps, cWps * Yps, xcenter = FALSE, ycenter = FALSE)
    nomy <- names(ccas$ycoef[, 1])
    nomx <- names(ccas$xcoef[, 1])
    nlx <- nrow(Xps[, colnames(Xps) %in% nomx])
    ncx <- ncol(Xps[, colnames(Xps) %in% nomx])
    Xps2 <- Xps[, colnames(Xps) %in% nomx]
    Yps <- NULL

    if (d == 1) {
      # Calculate Canonical variables
      echap <- (Xps2) %*% ccas$xcoef[, 1]
      Xps2 <- NULL
      if (!is.null(X2)) {
        Xpt <- t((t(X2) - mdX) / sdX)
        Xpt2 <- Xpt[, colnames(Xps) %in% nomx]
        Xpt <- NULL
        Xps <- NULL
        mdX <- NULL
        sdX <- NULL
        echap2 <- (Xpt2) %*% ccas$xcoef[, 1]
        Xpt2 <- NULL
        echaplist <- list(CV_X_A = echap, CV_X_B = echap2, CCA = ccas)
        return(echaplist)
      }
      else {
        echaplist <- list(CV_X_A = echap, Y_A_cv=cbind(Y,echap=echap), CVec_X_B = ccas$xcoef, CCA = ccas)
        return(echaplist)
      }
    }
    else if (d >= 2) {
      # Calculate Canonical variables
      #echap <- cbind((Xps2) %*% ccas$xcoef[, 1], (Xps2) %*% ccas$xcoef[, 2])
      echap <- (Xps2) %*% ccas$xcoef[, 1:d]
      Xps2 <- NULL
      if (!is.null(X2)) {
        Xpt <- t((t(X2) - mdX) / sdX)
        Xpt2 <- Xpt[, colnames(Xps) %in% nomx]
        Xpt <- NULL
        Xps <- NULL
        mdx <- NULL
        sdX <- NULL
        #echap2 <- cbind((Xpt2) %*% ccas$xcoef[, 1], (Xpt2) %*% ccas$xcoef[, 2])
        echap2 <- (Xpt2) %*% ccas$xcoef[, 1:d]
        Xpt2 <- NULL
        echaplist <- list(CV_X_A = echap, CV_X_B = echap2, CCA = ccas)
        return(echaplist)
      }
      else {
        echaplist <- list(CV_X_A = echap, Y_A_cv=cbind(Y,echap=echap), CVec_X_B = ccas$xcoef, CCA = ccas)
        return(echaplist)
      }
    }
    else {
      print("d must be a positive integer")
    }
  }
}
