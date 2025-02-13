#' @title statMatch.superOM.phase2
#' @description Statistical matching (phase2) using superOM.
#' @param df.rec a data.frame, the recipient data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param opts a list, supplied using the function \code{statMatch.superOM.options()}, containing all the options.
#' @return A list, ...
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @noRd

statMatch.superOM.phase2 <- function(df.rec, df.don, don.weights, names.CV, names.NCV, zero.constraints, opts) {

  if (opts$print.details) {
    message("\nPhase 2 - Matching continuous variables.")
    message("\nCommon variables: ", paste0(names.CV, collapse = ", "), ".")
    message("\nNon-common variables: ", paste0(names.NCV, collapse = ", "), ".")
  }

  # Initialize grid of parameters and perform tuning of the hyperparameters via k-fold cross-validation
  if (nrow(opts$p2_params) != 1)
    tune <- statMatch.superOM.phase2.tuning(df.don, names.CV, names.NCV, don.weights, zero.constraints, opts$p2_params, opts)
  else tune <- opts$p2_params

  if (opts$print.details)
    message("\nPhase 2. Fit and predict tuned model ...")

  # Prepare data
  data <- list()
  data$mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)))
  data$mtx.don.NCV.raw <- df2mtx(df.don %>% select(all_of(names.NCV)))
  data$mtx.rec.CV.raw <- df2mtx(df.rec %>% select(all_of(names.CV)))
  if (opts$scaling != "no")
    data <- statMatch.scale(data, don.weights, opts$scaling)
  # Fit and predict using the best hyperparameters
  set.seed(opts$seed_phase2)
  fit <- statMatch.superOM.fit(data, tune$xdim[1], tune$ydim[1], tune$weight.CV[1], don.weights, opts)
  pred_raw <- statMatch.superOM.phase2.predict(fit, data, don.weights, opts$scaling, zero.constraints, keep_unconstrained = T)
  pred <- mtx2df(pred_raw, var_names = names.NCV)

  out <- list(tune = tune,
              fit = fit,
              df.match.don = pred,
              df.match = cbind(df.rec, pred))
}
