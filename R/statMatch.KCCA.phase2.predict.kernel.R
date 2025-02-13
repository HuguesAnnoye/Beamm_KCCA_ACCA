#' @title Predict for KCCA (or CCA)
#' @description Do prediction.
#' @param uX a matrix with the canonical variables for data set 1.
#' @param uY a matrix with the canonical variables for data set 2.
#' @param weights a vector of individual weights.
#' @param kernel_predict Type of kernel use for prediction.
#' @param d a positive integer, number of latent variable used in CCA.
#' @param h a numerical value, the bandwidth.
#' @param rot a logical, if TRUE the bandwidth is h multiply by the variance
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param matrix.tot.possibilitities Matrix containing the compatibilities.
#' @param print.details a logical, if TRUE print the details.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom stats sd
#' @noRd

statMatch.KCCA.phase2.predict.kernel <- function(uX,
                                                 uY,
                                                 h = 1,
                                                 d = 1,
                                                 kernel_predict = c("gauss", "gaussmin", "unif", "epan", "dist", "alea"),
                                                 rot = FALSE,
                                                 matrix.tot.possibilitities = NULL,
                                                 names.NCV = NULL,
                                                 weights = NULL,
                                                 print.details = FALSE) {
  if (ncol(uX) != ncol(uY)) {
    stop("Error :: the canonical variable have not the same dimension")
  }
  # Rule of thumb
  if (rot == TRUE) {
    h2 <- h * sd(uX)
  } else {
    h2 <- h
  }
  # if there is no matrix of compatibilities
  if (is.null(matrix.tot.possibilitities)) {
    matrix.tot.possibilitities <- matrix(1L, nrow = nrow(uX), ncol = nrow(uY))
  }
  # Weights
  if (is.null(weights)) {
    Wps <- rep(1, nrow(uX))
  } else {
    Wps <- weights
  }
  # Gaussian kernel
  if (kernel_predict == "gauss") {
    # Dimension
    if (d == 1) {
      if (print.details) message("d = 1")
      hmin <- compute_hmin_kcca_d1_rcpp(uX, uY, h)
      luY <- length(uY)
      luX <- length(uX)
      W <- compute_omega_kcca_d1_rcpp(uX, uY, hmin)
      if (print.details) message("nrow(W) : ", nrow(W), " and ncol(W) : ", ncol(W))
      # Calculate Omega
      W2 <- matrix.tot.possibilitities * W * Wps
      W <- NULL
      Wps <- NULL
      # if d=2
    } else if (d >= 2) {
      if (print.details) message(paste0("d = ", d))
      hmin <- compute_hmin_kcca_rcpp(uX, uY, h)
      W <- compute_omega_kcca_rcpp(uX, uY, hmin)
      W2 <- matrix.tot.possibilitities * W * Wps
    } else {
      message("d must be positive")
    }
  } else if (kernel_predict == "alea") {
    # new methodologie
    if (d == 1) {
      if (print.details) message("d = 1")
      hmin <- compute_hmin_with_comp_mtx_kcca_d1_rcpp(uX, uY, matrix.tot.possibilitities, h)
      luY <- length(uY)
      luX <- length(uX)
      W <- compute_omega_kcca_d1_rcpp(uX, uY, hmin)
      if (print.details) message("nrow(W) : ", nrow(W), " and ncol(W) : ", ncol(W))
      matrix.tot <- matrix.tot.possibilitities
      W2 <- matrix.tot * W * Wps
      matrix.tot <- NULL
      W <- NULL
      Wps <- NULL
    } else if (d >= 2) {
      if (print.details) message(paste0("d = ", d))
      hmin <- compute_hmin_with_comp_mtx_kcca_rcpp(uX, uY, matrix.tot.possibilitities, h)
      W <- compute_omega_kcca_rcpp(uX, uY, hmin)
      if (print.details) message("nrow(W) : ", nrow(W), " and ncol(W) : ", ncol(W))
      W2 <- matrix.tot.possibilitities * W * Wps
      W <- NULL
      Wps <- NULL
    } else {
      message("d must be positive")
    }
  } else if (kernel_predict == "epan") {
    if (d == 1) {
      if (print.details) message("d = 1")
      hc <- expand.grid(X = uX, Y = uY)
      luY <- length(uY)
      luX <- length(uX)
      h3 <- compute_hmin_kcca_d1_rcpp(uX, uY, h)
      h4 <- rep(h3, each = luX) # à vérifier
      h3 <- NULL
      hc <- cbind(hc, h4)
      h4 <- NULL
      echap5 <- (hc[1] - hc[2]) / hc[3]
      if (print.details) message("echap5 : OK")
      hc <- NULL
      Kern <- apply(echap5, 1, function(u) epan_new(u)) # apply(hc, 1, function(u) gausskcca_new(x = u[1], y = u[2], h = u[3]))
      if (print.details) message("Kern : OK")
      echap5 <- NULL
      W <- matrix(Kern, nrow = luX, ncol = luY)
      Kern <- NULL
      if (print.details) message("nrow(W) : ", nrow(W), " and ncol(W) : ", ncol(W))
      W2 <- matrix.tot.possibilitities * W * Wps
      matrix.tot.possibilitities <- NULL
      W <- NULL
      Wps <- NULL
    } else if (d == 2) {
      if (print.details) message("d = 2")
      h3 <- compute_hmin_kcca_rcpp(uX, uY, h)
      h4 <- rep(h3, each = nrow(uX))
      h3 <- NULL
      echap3 <- expand.grid(X = uX[, 1], Y = uY[, 1])
      echap3 <- cbind(echap3, h4)
      echap4 <- expand.grid(Xdi = uX[, 2], Y = uY[, 2])
      echap4 <- cbind(echap4, h4)
      h4 <- NULL
      echap5 <- (echap3[, 1] - echap3[, 2]) / echap3[, 3]
      echap6 <- (echap4[, 1] - echap4[, 2]) / echap3[, 3]
      echap3 <- NULL
      echap4 <- NULL
      hc <- cbind(echap5, echap6)
      echap5 <- NULL
      echap6 <- NULL
      Kern <- apply(hc, 1, function(u) epan3d(u))
      if (print.details) message("Kern : OK")
      hc <- NULL
      W <- matrix(Kern, nrow = nrow(uX), ncol = nrow(uY))
      uX <- NULL
      uY <- NULL
      if (print.details) message("nrow(W) : ", nrow(W), " and ncol(W) : ", ncol(W))
      W2 <- matrix.tot.possibilitities * W * Wps
      W <- NULL
      Wps <- NULL
    } else {
      message("d must be <  3 and > 0")
    }
  }
  return(W2)
}
