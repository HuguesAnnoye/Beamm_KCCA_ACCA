#' @title statMatch.CCA.phase2
#' @description Statistical matching (phase2) using CCA
#' @param df.rec a data.frame, the receiver data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param opts a list returned by the function \code{statMatch.CCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#'
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr '%>%'
#' @noRd

statMatch.CCA.phase2 <- function(df.rec, df.don, don.weights, names.CV, names.NCV, zero.constraints, opts) {

  if (opts$print.details) {
    message("\nPhase 2 - Matching continuous variables.")
    message("\nCommon variables: ", paste0(names.CV, collapse = ", "), ".")
    message("\nNon-common variables: ", paste0(names.NCV, collapse = ", "), ".")
  }

  # Tune the hyperparameters (if necessary)
  if (opts$P2$n_h != 1) {
    if (isTRUE(opts$par)) {
      tune <- statMatch.CCA.phase2.tuning.par(df.don, names.CV, names.NCV, don.weights, zero.constraints, opts)
    } else {
      tune <- statMatch.CCA.phase2.tuning(df.don, names.CV, names.NCV, don.weights, zero.constraints, opts)
    }
  } else {  tune <- list(results = opts$P2$hpar)}

  if (opts$print.details)
    message("\nPhase 2. Fit and predict tuned model ...")

  # Prepare data
  data <- list()
  data$mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)), n_1_levels = TRUE)
  data$mtx.don.NCV.raw <- df2mtx(df.don %>% select(all_of(names.NCV)), n_1_levels = TRUE)
  data$mtx.rec.CV.raw <- df2mtx(df.rec %>% select(all_of(names.CV)), n_1_levels = TRUE)
  data$mtx.don.NCV.n <- df2mtx(df.don %>% select(all_of(names.NCV)), n_1_levels = FALSE)
  #if (opts$scaling != "no") No scaling for CCA because the scaling is mendatory
    #data <- statMatch.scale(data, don.weights, opts$scaling)

  # Fit and predict using the best hyperparameters

  fit <- statMatch.CCA.fit(data$mtx.don.CV.raw,
                            data$mtx.don.NCV.raw,
                            weights = don.weights,
                            X2 = data$mtx.rec.CV.raw,
                            d = tune$results$d[1])

  # Compute the predictions
  if (opts$type_predict == "matrix") {
    # Do the prediction using a matrix
  pred_raw <- statMatch.KCCA.phase2.predict(CV_X_A = fit$CV_X_A, #canonical variables for data set A
                                            CV_X_B = fit$CV_X_B, #canonical variables for data set B
                                            Y = data$mtx.don.NCV.n, #Non-common variables for data set A
                                            h = tune$results$h[1],
                                            d = tune$results$d[1],
                                            kernel_predict = opts$P2$kernel_predict,
                                            rot = opts$rot,
                                            names.NCV = colnames(data$mtx.don.NCV.n),
                                            weights = don.weights,
                                            data_cat_rec = df.rec,
                                            data_cat_don = df.don,
                                            zero.constraints = zero.constraints,
                                            matrix.tot.possibilitities = NULL,
                                            print.details = opts$print.details)
  } else if (opts$type_predict == "loop") {
    # Do the prediction using a rcpp loop
    data$mtx.don.CV.n <- df2mtx(df.don %>% select(all_of(names.CV)), n_1_levels = FALSE)
    data$mtx.rec.CV.n <- df2mtx(df.rec %>% select(all_of(names.CV)), n_1_levels = FALSE)
    pred_raw <- statMatch.CCA.phase2.predict(CV_X_A = fit$CV_X_A, #canonical variables for data set A
                                              CV_X_B = fit$CV_X_B, #canonical variables for data set B
                                              h = tune$results$h[1],
                                              d = tune$results$d[1],
                                              mtx.rec.CV.raw = data$mtx.rec.CV.n,
                                              mtx.don.CV.raw = data$mtx.don.CV.n,
                                              mtx.don.NCV.raw = data$mtx.don.NCV.n,
                                              don.weights = don.weights,
                                              opts = opts,
                                              zero.constraints = zero.constraints,
                                              keep_unconstrained = T)
    colnames(pred_raw) <- colnames(data$mtx.don.NCV.n)
  }
  pred <- mtx2df(pred_raw, var_names = names.NCV)
  out <- list(tune = tune,
              fit = fit,
              df.match.don = pred,
              df.match = cbind(df.rec, pred))
  return(out)
}
