#' @title statMatch.PCA.phase2
#' @description Statistical matching (phase2) using PCA-CCA
#' @param df.rec a data.frame, the receiver data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param opts a list returned by the function \code{statMatch.PCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#'
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr '%>%'
#' @noRd

statMatch.PCA.phase2 <- function(df.rec, df.don, don.weights, names.CV, names.NCV, zero.constraints, opts) {
  if (opts$print.details) {
    message("\nPhase 2 - Matching continuous variables.")
    message("\nCommon variables: ", paste0(names.CV, collapse = ", "), ".")
    message("\nNon-common variables: ", paste0(names.NCV, collapse = ", "), ".")
  }

  # Tune the hyperparameters (if necessary)
  if ((opts$P2$n_combs) * (opts$P2$n_h) != 1) {
    tune <- statMatch.PCA.phase2.tuning(df.don, names.CV, names.NCV, don.weights, zero.constraints, opts)
  } else {
    tune <- list(results = opts$P2$hpar)
  }

  if (opts$print.details) {
    message("\nPhase 2. Fit and predict tuned model ...")
  }

  # Prepare data
  data <- list()
  data$mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)))
  data$mtx.don.NCV.raw <- df2mtx(df.don %>% select(all_of(names.NCV)))
  data$mtx.rec.CV.raw <- df2mtx(df.rec %>% select(all_of(names.CV)))
  if (opts$scaling != "no") {
    data <- statMatch.scale(data, don.weights, opts$scaling)
  }

  # Fit and predict using the best hyperparameters

  fit <- statMatch.PCA.fit(data$mtx.don.CV.raw,
    data$mtx.don.NCV.raw,
    weights = don.weights,
    X2 = data$mtx.rec.CV.raw,
    d = tune$results$d[1],
    lat_X = tune$results$lat_X[1],
    lat_Y = tune$results$lat_Y[1]
  )

  # Compute the predictions
  if (opts$type_predict == "matrix") {
    # Do the prediction using a matrix
    pred_raw <- statMatch.KCCA.phase2.predict(
      CV_X_A = fit$CV_X_A, # canonical variables for data set A
      CV_X_B = fit$CV_X_B, # canonical variables for data set B
      Y = data$mtx.don.NCV.raw, # Non-common variables for data set A
      h = tune$results$h[1],
      d = tune$results$d[1],
      kernel_predict = opts$P2$kernel_predict,
      scaling = opts$scaling,
      rot = opts$rot,
      names.NCV = colnames(data$mtx.don.NCV.raw),
      weights = don.weights,
      data_cat_rec = df.rec,
      data_cat_don = df.don,
      zero.constraints = zero.constraints,
      matrix.tot.possibilitities = NULL,
      print.details = opts$print.details
    )
  } else if (opts$type_predict == "loop") {
    # Do the prediction using a rcpp loop
    pred_raw <- statMatch.CCA.phase2.predict(
      CV_X_A = fit$CV_X_A, # canonical variables for data set A
      CV_X_B = fit$CV_X_B, # canonical variables for data set B
      h = tune$results$h[1],
      d = tune$results$d[1],
      mtx.rec.CV.raw = data$mtx.rec.CV.raw,
      mtx.don.CV.raw = data$mtx.don.CV.raw,
      mtx.don.NCV.raw = data$mtx.don.NCV.raw,
      don.weights = don.weights,
      opts = opts,
      zero.constraints = zero.constraints,
      keep_unconstrained = T
    )
    colnames(pred_raw) <- colnames(data$mtx.don.NCV.raw)
  }
  pred <- mtx2df(pred_raw, var_names = names.NCV)
  out <- list(
    tune = tune,
    fit = fit,
    df.match.don = pred,
    df.match = cbind(df.rec, pred)
  )
  return(out)
}
