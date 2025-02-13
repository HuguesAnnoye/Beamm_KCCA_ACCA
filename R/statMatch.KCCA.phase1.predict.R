#' @title Predict for CCA and KCCA
#' @description Do prediction
#' @param CV_X_A a matrix with the canonical variables for data set A.
#' @param CV_X_B a matrix with the canonical variables for data set B.
#' @param Y a matrix with the non-common variables in data set A.
#' @param weights a vector of individual weights.
#' @param kernel_predict a character string, the type of kernel use for prediction.
#' @param d a positive integer, number of latent variable used in KCCA.
#' @param h a numerical value, Bandwidth.
#' @param rot a logical, if TRUE the bandwidth is h multiply by the variance.
#' @param matrix.tot.possibilitities a numeric matrix containing the compatibilities.
#' @param print.details a logical, if TRUE print the details.
#' @param names.NCV a character vector, with the names of the non-common variables used for statistical matching.
#' @param scaling a character string, the name of the method used to scale the variables in the receiver and donor data sets: \code{"z-score"}, \code{"min-max"} or \code{"no"} (no scaling).
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom stats sd
#' @noRd


statMatch.KCCA.phase1.predict <- function(CV_X_A,
                                          CV_X_B,
                                          Y,
                                          h = 1,
                                          d = 1,
                                          kernel_predict = c("gauss", "gaussmin", "unif", "epan", "dist", "alea"),
                                          scaling = c("no", "z-score", "min-max"),
                                          # X1 = NULL,
                                          # X2 = NULL,
                                          rot = FALSE,
                                          matrix.tot.possibilitities = NULL,
                                          names.NCV = NULL,
                                          weights = NULL,
                                          print.details = FALSE) {
  if (ncol(CV_X_A) != ncol(CV_X_B)) {
    stop("Error :: the canonical variable have not the same dimension")
  }

  if (rot == TRUE) {
    h2 <- h * sd(CV_X_A)
  }
  else {
    h2 <- h
  }
  if (is.null(names.NCV)) {
    names.NCV <- colnames(Y)
  }

  if (is.null(matrix.tot.possibilitities)) {
    matrix.tot.possibilitities <- matrix(1, nrow = nrow(CV_X_A), ncol = nrow(CV_X_B))
    # Wps <- rep(1,nrow(CV_X_A))
  }

  if (is.null(weights)) {
    Wps <- rep(1, nrow(CV_X_A))
  } else {
    Wps <- weights
  }

  if (kernel_predict == "gauss") {
    if (d == 1) {
      if (print.details) message("d = 1")
      luY <- length(CV_X_B)
      luX <- length(CV_X_A)
      hmin <- compute_hmin_with_comp_mtx_kcca_d1_rcpp(CV_X_A, CV_X_B, matrix.tot.possibilitities, h)
      W <- compute_omega_kcca_d1_rcpp(CV_X_A, CV_X_B, hmin)
      if (print.details) message("nrow(W) : ", nrow(W), " and ncol(W) : ", ncol(W))
      W2 <- matrix.tot.possibilitities * W * Wps
      # matrix.tot <- NULL
      W <- NULL
      Wps <- NULL
      Wcs <- colSums(W2)
      W1 <- t(W2) / Wcs
      W2 <- NULL
      Wcs <- NULL
      data_testchap <- W1 %*% Y
      W1 <- NULL
      mtx_PRED <- data_testchap[, names.NCV, drop = F]
    } else if (d >= 2) {
      if (print.details) message(paste0("d = ", d))
      hmin <- compute_hmin_with_comp_mtx_kcca_rcpp(CV_X_A, CV_X_B, matrix.tot.possibilitities, h)
      W <- compute_omega_kcca_rcpp(CV_X_A, CV_X_B, hmin)
      W2 <- matrix.tot.possibilitities * W * Wps
      W <- NULL
      Wps <- NULL
      Wcs <- colSums(W2)
      W1 <- t(W2) / Wcs
      W2 <- NULL
      Wcs <- NULL
      data_testchap <- W1 %*% Y
      W1 <- NULL
      mtx_PRED <- data_testchap[, names.NCV, drop = F]
    } else {
      message("d must be >  0")
    }
  } else if (kernel_predict == "alea") {
    # new methodologie
    if (d == 1) {
      if (print.details) message("d = 1")
      hmin <- compute_hmin_with_comp_mtx_kcca_d1_rcpp(CV_X_A, CV_X_B, matrix.tot.possibilitities, h)
      luY <- length(CV_X_B)
      luX <- length(CV_X_A)
      W <- compute_omega_kcca_d1_rcpp(CV_X_A, CV_X_B, hmin)
      if (print.details) message("nrow(W) : ", nrow(W), " and ncol(W) : ", ncol(W))
      matrix.tot <- matrix.tot.possibilitities
      # W2 <- matrix.tot[,] * W * Wps
      W2 <- matrix.tot * W * Wps
      matrix.tot <- NULL
      W <- NULL
      Wps <- NULL
      Wcs <- colSums(W2)
      W1 <- t(W2) / Wcs
      W2 <- NULL
      Wcs <- NULL
      data_testchap <- matrix(NA, nrow = luY, ncol = length(names.NCV))
      colnames(data_testchap) <- names.NCV
      for (j in 1:luY) {
        iraw <- sample(1:luX, 1, prob = W1[j, ])
        # message("iraw = ",iraw)
        # message(names.NCV)
        data_testchap[j, ] <- Y[iraw, names.NCV, drop = F]
      }
      W1 <- NULL
      mtx_PRED <- data_testchap
    } else if (d >= 2) {
      if (print.details) message(paste0("d = ", d))
      hmin <- compute_hmin_with_comp_mtx_kcca_rcpp(CV_X_A, CV_X_B, matrix.tot.possibilitities, h)
      W <- compute_omega_kcca_rcpp(CV_X_A, CV_X_B, hmin)
      W2 <- matrix.tot.possibilitities * W * Wps
      W <- NULL
      Wps <- NULL
      Wcs <- colSums(W2)
      W1 <- t(W2) / Wcs
      W2 <- NULL
      Wcs <- NULL
      luY <- nrow(CV_X_B)
      luX <- nrow(CV_X_A)
      data_testchap <- matrix(NA, nrow = luY, ncol = length(names.NCV))
      colnames(data_testchap) <- names.NCV
      for (j in 1:luY) {
        iraw <- sample(1:luX, 1, prob = W1[j, ])
        data_testchap[j, ] <- Y[iraw, names.NCV, drop = F]
      }
      W1 <- NULL
      mtx_PRED <- data_testchap
    } else {
      message("d must be > 0")
    }
  } else if (kernel_predict == "epan") {
    if (d == 1) {
      if (print.details) message("d = 1")
      hc <- expand.grid(X = CV_X_A, Y = CV_X_B)
      luY <- length(CV_X_B)
      luX <- length(CV_X_A)
      h3 <- compute_hmin_with_comp_mtx_kcca_d1_rcpp(CV_X_A, CV_X_B, matrix.tot.possibilitities, h)
      h4 <- rep(h3, each = luX) # à vérifier
      # uX <- NULL
      h3 <- NULL
      hc <- cbind(hc, h4)
      h4 <- NULL
      # echap5 <- as.big.matrix((hc[1] - hc[2]) / hc[3])
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
      Wcs <- colSums(W2)
      W1 <- t(W2) / Wcs
      W2 <- NULL
      Wcs <- NULL
      data_testchap <- W1 %*% Y
      W1 <- NULL
      mtx_PRED <- data_testchap[, names.NCV, drop = F]
    } else if (d == 2) {
      if (print.details) message("d = 2")
      h3 <- compute_hmin_with_comp_mtx_kcca_rcpp(CV_X_A, CV_X_B, matrix.tot.possibilitities, h)
      hmin2 <- NULL
      h4 <- rep(h3, each = nrow(CV_X_A)) # à vérifier
      h3 <- NULL
      echap3 <- expand.grid(X = CV_X_A[, 1], Y = CV_X_B[, 1])
      echap3 <- cbind(echap3, h4)
      echap4 <- expand.grid(Xdi = CV_X_A[, 2], Y = CV_X_B[, 2])
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
      W <- matrix(Kern, nrow = nrow(CV_X_A), ncol = nrow(CV_X_B))
      CV_X_A <- NULL
      CV_X_B <- NULL
      if (print.details) message("nrow(W) : ", nrow(W), " and ncol(W) : ", ncol(W))
      W2 <- matrix.tot.possibilitities * W * Wps
      W <- NULL
      Wps <- NULL
      Wcs <- colSums(W2)
      W1 <- t(W2) / Wcs
      W2 <- NULL
      Wcs <- NULL
      data_testchap <- W1 %*% Y
      W1 <- NULL
      mtx_PRED <- data_testchap[, names.NCV, drop = F]
    } else {
      message("d must be <  2")
    }
  }

  # Go back to the original scale of the variables
  if (scaling == "z-score") {
    mtx_PRED <- scale_mtx_inv(mtx_PRED,
                              center = attr(Y, "scaled:center"),
                              scale = attr(Y, "scaled:scale"))
  } else if (scaling == "min-max") {
    mtx_PRED <- normalize_mtx_inv(mtx_PRED,
                                  min = attr(Y, "min"),
                                  max = attr(Y, "max"))
  }

  return(mtx_PRED)
}
