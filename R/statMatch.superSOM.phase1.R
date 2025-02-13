#' @title statMatch.superOM.phase1
#' @description Statistical matching (phase1) using superOM.
#' @param df.rec a data.frame, the recipient data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param opts a list, supplied using the function \code{statMatch.superOM.options()}, containing all the options.
#' @return A list, ...
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @noRd

statMatch.superOM.phase1 <- function(df.rec, df.don, don.weights, names.CV, names.NCV, opts) {

  # Initialization
  names.CV.cat <- names.CV[unlist(lapply(df.don %>% select(all_of(names.CV)), is.factor))]
  names.NCV.cat <- names.NCV[unlist(lapply(df.don %>% select(all_of(names.NCV)), is.factor))]

  if (length(names.CV.cat) > 0 & length(names.NCV.cat) > 0) {
    if (opts$print.details) {
      message("\nPhase 1 - Matching categorical variables.")
      message("\nCommon variables: ", paste0(names.CV, collapse = ", "), ".")
      message("\nNon-common variables: ", paste0(names.NCV.cat, collapse = ", "), ".")
    }

    # Initialize grid of parameters and perform tuning of the hyperparameters via k-fold cross-validation
    if (nrow(opts$p1_params) != 1)
      tune <- statMatch.superOM.phase1.tuning(df.don, names.CV, names.NCV.cat, don.weights, opts$p1_params, opts)
    else tune <- opts$p1_params

    if (opts$print.details)
      message("\nPhase 1. Fit and predict tuned model ...")

    # Prepare data
    data <- list()
    data$mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)))
    data$mtx.don.NCV.raw <- df2mtx(df.don %>% select(all_of(names.NCV.cat)))
    data$mtx.rec.CV.raw <- df2mtx(df.rec %>% select(all_of(names.CV)))
    if (opts$scaling != "no")
      data <- statMatch.scale(data, don.weights, opts$scaling)

    # Fit and predict using the best hyperparameters
    set.seed(opts$seed_phase1)
    fit <- statMatch.superOM.fit(data, tune$xdim[1], tune$ydim[1], tune$weight.CV[1], don.weights, opts)
    comp.mtx <- compute_mtx_compatibilities(df.rec, df.don, names.CV.cat, opts$print.details, index = F)
    pred_raw <- statMatch.superOM.phase1.predict(fit, data, don.weights, comp.mtx, opts$scaling)
    pred <- mtx2df(pred_raw, var_names = names.NCV.cat)

    # Compute the sum of distances between the predictions and their compatible vectors in the donor data set
    #   This vector is used when averaging multiple models, for stability reasons.
    sum_dist_comp_donor <- compute_sum_dist_comp_donor(pred_raw, data$mtx.don.NCV.raw, comp.mtx)

    out <- list(tune = tune,
                fit = fit,
                df.match.don = pred,
                df.match = cbind(df.rec, pred),
                sum_dist_comp_donor = sum_dist_comp_donor)
  } else {
    if (length(names.NCV.cat) == 0)
      message("\nPhase 1 - No categorical variables to predict. Passing to the next phase ...")
    if (length(names.CV.cat) == 0 & length(names.NCV.cat) > 0)
      stop("No common categorical variables. Impossible to compute compatibility matrix.")
    out <- list(df.match = df.rec)
  }
  return(out)
}


compute_sum_dist_comp_donor <- function(pred, data_don, comp_mtx) {
  n_rec <- nrow(pred)
  out <- rep(NA, n_rec)
  for (i in seq_len(n_rec)) {
    out[i] <- sum(compute_eucl_dist_vec_mat_rcpp(pred[i,], data_don[comp_mtx[,i] == 1L, , drop = F]))
  }
  return(out)
}
