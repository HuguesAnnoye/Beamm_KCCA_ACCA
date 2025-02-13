#' @title Fitting for KCCA
#' @description Fit with KCCA.
#' @param X a matrix with the common variables in data set 1.
#' @param Y a matrix with the non-common variables in data set 1.
#' @param X2 a matrix with the common variables in data set 2.
#' @param weights a vector of individual weights.
#' @param d a positive integer, number of latent variable used in CCA.
#' @param h Value for h, the bandwidth of the KCCA kernel.
#' @param hy Second value for h, if we want data tha bandwidth will be different for the non-common variable.
#' @param g Value for g, the regularization parameter in KCCA.
#' @param rot a logical, if TRUE the bandwidths are multiply by the variance.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @rawNamespace import(kernlab, except = predict)
#' @importFrom methods is
#' @importFrom stats sd
#' @noRd


statMatch.KCCA.fit <- function(X,
                               Y,
                               X2 = NULL,
                               weights = NULL,
                               d = 1,
                               h = 1,
                               hy = NULL,
                               g = 1,
                               rot = FALSE) {
  if (is.null(weights)) {

    ## Multiply the bandwidth by the variance if rot == TRUE
    if (rot == TRUE) {
      pb <- function(x, y) {
        (x - y) %*% (x - y)
      }
      mdist <- function(x, y) {
        t(apply(x, 1, function(p) apply(y, 1, function(q) pb(p, q))))
      }
      hs <- h / sd(mdist(X, X))
      if (is.null(hy)) {
        hys <- hs
      }
      else {
        hys <- hy / sd(mdist(Y, Y))
      }
      # KCCA
      kccas <- kccaw(X, Y, kpar = list(sigma = hs), kpar2 = list(sigma = hys), gamma = g, ncomps = 2 * d - 1)
      kpar <- list(sigma = hs)
      kpar2 <- list(sigma = hys)
    } else {
      if (is.null(hy)) {
        hy <- h
      }
      ## KCCA
      kccas <- kccaw(X, Y, kpar = list(sigma = h), kpar2 = list(sigma = hy), gamma = g, ncomps = 2 * d - 1)
      kpar <- list(sigma = h)
      kpar2 <- list(sigma = hy)
    }
    # create the kernel for the prediction part
    kernel <- "rbfdot"
    if (!is(kernel, "kernel")) {
      if (is(kernel, "function")) {
        kernel1 <- deparse(substitute(kernel))
        # kernel2 <- deparse(substitute(kernel))
      }
      kernel1 <- do.call(kernel, kpar)
    }
    if (!is(kernel1, "kernel")) stop("kernel must inherit from class `kernel'")
    kxx <- kernelMatrix(kernel1, x = as.matrix(X))
    if (d == 1) {
      # Calculate the canonical variables for the donnor
      echapk <- (kxx) %*% kccas@xcoef[, 1]
      if (!is.null(X2)) {
        # Calculate the canonical variables for the receiver
        kxy <- kernelMatrix(kernel1, x = as.matrix(X2), y = as.matrix(X))
        echapk2 <- (kxy) %*% kccas@xcoef[, 1]
        echaplist <- list(CV_X_A = echapk, CV_X_B = echapk2, CCA = kccas) # list(echapk, echapk2)
        return(echaplist)
      }
      else {
        echaplist <- list(CV_X_A = echapk, CVec_X_B = kccas@xcoef, CCA = kccas)
        return(echaplist)
      }
    } else if (d >= 2) {
      # Calculate the canonical variables
      # echapk <- cbind((kxx) %*% kccas@xcoef[, 1], (kxx) %*% kccas@xcoef[, 3])
      echapk <- (kxx) %*% kccas@xcoef[, seq(1, 2 * d - 1, 2)]
      if (!is.null(X2)) {
        # Calculate the canonical variables for the receiver
        kxy <- kernelMatrix(kernel1, x = as.matrix(X2), y = as.matrix(X))
        # echapk2 <- cbind((kxy) %*% kccas@xcoef[, 1], (kxy) %*% kccas@xcoef[, 3])
        echapk2 <- (kxy) %*% kccas@xcoef[, seq(1, 2 * d - 1, 2)]
        echaplist <- list(CV_X_A = echapk, CV_X_B = echapk2, CCA = kccas) # list(echapk, echapk2)
        return(echaplist)
      }
      else {
        echaplist <- list(CV_X_A = echapk, CVec_X_B = kccas@xcoef, CCA = kccas)
        return(echaplist)
      }
    } else {
      print("d must be a positive integer")
    }
  } else {
    Wps <- weights

    ## KCCA
    if (rot == TRUE) {
      pb <- function(x, y) {
        (x - y) %*% (x - y)
      }
      mdist <- function(x, y) {
        t(apply(x, 1, function(p) apply(y, 1, function(q) pb(p, q))))
      }

      ## Multiply the bandwidth by the variance if rot = TRUE
      hs <- h / sd(mdist(X, X))
      if (is.null(hy)) {
        hys <- hs
      }
      else {
        hys <- hy / sd(mdist(Y, Y))
      }
      ## KCCA
      kccas <- kccaw(X, Y, w = Wps, kpar = list(sigma = hs), kpar2 = list(sigma = hys), gamma = g, ncomps = 2 * d - 1)
      kpar <- list(sigma = hs)
      kpar2 <- list(sigma = hys)
    } else {
      if (is.null(hy)) {
        hy <- h
      }
      kccas <- kccaw(X, Y, w = Wps, kpar = list(sigma = h), kpar2 = list(sigma = hy), gamma = g, ncomps = 2 * d - 1)
      kpar <- list(sigma = h)
      kpar2 <- list(sigma = hy)
    }
    # create the kernel for the prediction part
    kernel <- "rbfdot"
    if (!is(kernel, "kernel")) {
      if (is(kernel, "function")) {
        kernel1 <- deparse(substitute(kernel))
      }
      kernel1 <- do.call(kernel, kpar)
    }
    if (!is(kernel1, "kernel")) stop("kernel must inherit from class `kernel'")
    kxx <- kernelMatrix(kernel1, x = as.matrix(X))
    if (d == 1) {
      # Calculate the canonical variables for the donnor
      echapk <- (kxx) %*% kccas@xcoef[, 1]
      if (!is.null(X2)) {
        # Calculate the canonical variables for the receiver
        kxy <- kernelMatrix(kernel1, x = as.matrix(X2), y = as.matrix(X))
        echapk2 <- (kxy) %*% kccas@xcoef[, 1]
        echaplist <- list(CV_X_A = echapk, CV_X_B = echapk2, CCA = kccas) # list(echapk, echapk2)
        return(echaplist)
      } else {
        echaplist <- list(CV_X_A = echapk, CVec_X_B = kccas@xcoef, CCA = kccas)
        return(echaplist)
      }
    } else if (d >= 2) {
      # Calculate the canonical variables for the donnor
      # echapk <- cbind((kxx) %*% (kccas@xcoef[, 1]), (kxx) %*% (kccas@xcoef[, 3]))
      echapk <- cbind((kxx) %*% kccas@xcoef[, seq(1, 2 * d - 1, 2)])
      if (!is.null(X2)) {
        # Calculate the canonical variables for the receiver
        kxy <- kernelMatrix(kernel1, x = as.matrix(X2), y = as.matrix(X))
        # echapk2 <- cbind((kxy) %*% (kccas@xcoef[, 1]), (kxy) %*% (kccas@xcoef[, 3]))
        echapk2 <- cbind((kxy) %*% kccas@xcoef[, seq(1, 2 * d - 1, 2)])
        echaplist <- list(CV_X_A = echapk, CV_X_B = echapk2, CCA = kccas) # list(echapk, echapk2)
        return(echaplist)
      } else {
        echaplist <- list(CV_X_A = echapk, CVec_X_B = kccas@xcoef, CCA = kccas)
        return(echaplist)
      }
    } else {
      print("d must be a positive integer.")
    }
  }
}
