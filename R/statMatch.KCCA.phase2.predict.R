#' @title Predict for CCA and KCCA
#' @description Do prediction
#' @param CV_X_A a matrix, with the canonical variables for data set 1.
#' @param CV_X_B a matrix, with the canonical variables for data set 2.
#' @param Y a matrix, with the non-common variables in data set 1.
#' @param data_cat_rec a data frame, with the categorical variable for data set 1.
#' @param data_cat_don a data frame, with the categorical variable for data set 2.
#' @param zero.constraints Vector with the name of the zero constraints.
#' @param weights a vector of individual weights.
#' @param kernel_predict a character string, the type of kernel use for prediction
#' @param d a positive integer, number of latent variable used in KCCA.
#' @param h a numerical value, the bandwidth.
#' @param rot a logical, if TRUE the bandwidth is h multiply by the variance.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param matrix.tot.possibilitities a matrix, containing the compatibilities.
#' @param print.details a logical, if TRUE print the details.
#' @param scaling a character string, the name of the method used to scale the variables in the receiver and donor data sets: \code{"z-score"}, \code{"min-max"} or \code{"no"} (no scaling).
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom dplyr setdiff
#' @importFrom stats sd
#' @noRd

# attention il faut tenir compte du lenght variable =1
statMatch.KCCA.phase2.predict <- function(CV_X_A,
                                          CV_X_B,
                                          Y,
                                          data_cat_rec,
                                          data_cat_don,
                                          zero.constraints = NULL,
                                          h = 1,
                                          d = 1,
                                          kernel_predict = c("gauss", "unif", "epan", "dist", "alea"),
                                          scaling = c("no", "z-score", "min-max"),
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
  }

  if (is.null(weights)) {
    Wps <- rep(1, nrow(CV_X_A))
  } else {
    Wps <- weights
  }

  Wkernel <- statMatch.KCCA.phase2.predict.kernel(
    uX = CV_X_A,
    uY = CV_X_B,
    h = h,
    d = d,
    kernel_predict = kernel_predict,
    rot = rot,
    matrix.tot.possibilitities = matrix.tot.possibilitities,
    names.NCV = names.NCV,
    weights = weights,
    print.details = print.details
  )

  # Go back to the original scale of the variables
  if (scaling == "z-score") {
    Y <- scale_mtx_inv(Y,
      center = attr(Y, "scaled:center"),
      scale = attr(Y, "scaled:scale")
    )
  } else if (scaling == "min-max") {
    Y <- normalize_mtx_inv(Y,
      min = attr(Y, "min"),
      max = attr(Y, "max")
    )
  }


  names.NCV_without_ZC <- setdiff(names.NCV, zero.constraints)
  if (identical(names.NCV_without_ZC, character(0))) {
    mtx_PRED <- NULL
    name_PRED <- NULL
  } else {
    mtx_PRED <- statMatch.KCCA.phase2.predict.prediction(
      W_keep = Wkernel,
      Z = Y,
      d = d,
      kernel_predict = kernel_predict,
      matrix.tot.possibilitities = NULL,
      names.NCV = names.NCV_without_ZC,
      weights = weights,
      print.details = print.details
    )
    name_PRED <- colnames(mtx_PRED)
  }
  for (var in zero.constraints)
  {
    name_ZC <- paste0("ZC_", var)
    if (print.details) {
      message(name_ZC)
    }
    start.timezc <- Sys.time()
    matrix.tot.possibilitities <- compute_mtx_compatibilities(
      df.rec = data_cat_rec[, name_ZC, drop = F],
      df.don = data_cat_don[, name_ZC, drop = F],
      names.CV.cat = name_ZC,
      print.info = FALSE,
      index = FALSE
    )
    # matrix.tot.possibilitities <- statMatch.CCA.possibilities(data_cat_rec=data_cat_rec[ ,name_ZC,drop=F],
    #                                                                  data_cat_don=data_cat_don[ ,name_ZC,drop=F],
    #                                                                  names.CV=name_ZC,
    #                                                                  names.CV_Cat_prio=NULL)

    start.timezc2 <- Sys.time()
    if (print.details) {
      message(start.timezc)
      message(start.timezc2)
      message(start.timezc2 - start.timezc)
    }
    mtx_tmp <- statMatch.KCCA.phase2.predict.prediction(
      W_keep = Wkernel,
      Z = Y[, var, drop = F],
      d = d,
      kernel_predict = kernel_predict,
      matrix.tot.possibilitities = matrix.tot.possibilitities,
      weights = weights,
      names.NCV = var,
      print.details = print.details
    )

    mtx_PRED <- cbind(mtx_PRED, mtx_tmp)
    name_PRED <- c(name_PRED, var)
    colnames(mtx_PRED) <- name_PRED
    mtx_tmp <- NULL
  }

  return(mtx_PRED)
}
