#' @title CCA - PHASE 1 TUNING AND PREDICTION
#' @description Tuning and prediction for CCA statistical matching process
#' @param df.rec a data.frame, the receiver data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param opts a list returned by the function \code{statMatch.CCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr '%>%'
#' @importFrom keras k_clear_session
#' @importFrom tensorflow set_random_seed
#' @noRd


statMatch.CCA.phase1 <- function(df.rec, df.don, don.weights, names.CV, names.NCV, opts) {

  # Initialization
  names.CV.cat <- names.CV[unlist(lapply(df.don %>% select(all_of(names.CV)), is.factor))]
  names.NCV.cat <- names.NCV[unlist(lapply(df.don %>% select(all_of(names.NCV)), is.factor))]

  if (length(names.CV.cat) > 0 & length(names.NCV.cat) > 0) {
    if (opts$print.details) {
      message("\nPhase 1 - Matching categorical variables.")
      message("\nCommon variables: ", paste0(names.CV, collapse = ", "), ".")
      message("\nNon-common variables: ", paste0(names.NCV.cat, collapse = ", "), ".")
    }

    # Tune the hyperparameters (if necessary)
    if (opts$P1$n_h != 1) {
        if (isTRUE(opts$par)) {
          tune <- statMatch.CCA.phase1.tuning.par(df.don, names.CV, names.NCV.cat, don.weights, opts)
        } else {
          tune <- statMatch.CCA.phase1.tuning(df.don, names.CV, names.NCV.cat, don.weights, opts)
        }
    } else {
      tune <- list(results = opts$P1$hpar)
    }

    if (opts$print.details) {
      message("\nPhase 1. Fit and predict tuned model ...")
    }

    # Prepare data
    data <- list()
    data$mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)), n_1_levels = TRUE)
    data$mtx.don.NCV.raw <- df2mtx(df.don %>% select(all_of(names.NCV.cat)), n_1_levels = TRUE)
    data$mtx.rec.CV.raw <- df2mtx(df.rec %>% select(all_of(names.CV)), n_1_levels = TRUE)
    data$mtx.don.NCV.n <- df2mtx(df.don %>% select(all_of(names.NCV.cat)), n_1_levels = FALSE)
    # if (opts$scaling != "no") #The scaling is done in the fitting function
    #  data <- statMatch.scale(data, don.weights, opts$scaling)

    # Compute compatibility matrix
    comp.mtx <- compute_mtx_compatibilities(df.rec, df.don, names.CV.cat, opts$print.details, index = F)

    if (opts$print.details) message("Fitting")
    fit <- statMatch.CCA.fit(data$mtx.don.CV.raw,
      data$mtx.don.NCV.raw,
      weights = don.weights,
      X2 = data$mtx.rec.CV.raw,
      d = tune$results$d[1]
    )

    if (opts$print.details) message("Prediction")

    if (opts$type_predict == "matrix") {
      # Do the prediction using a matrix
      pred_raw <- statMatch.KCCA.phase1.predict(
        CV_X_A = fit$CV_X_A, # canonical variable for data set A
        CV_X_B = fit$CV_X_B, # canonical variable for data set B
        Y = data$mtx.don.NCV.n,
        h = tune$results$h[1],
        d = tune$results$d[1],
        kernel_predict = opts$P1$kernel_predict,
        rot = opts$rot,
        matrix.tot.possibilitities = comp.mtx,
        names.NCV = colnames(data$mtx.don.NCV.n),
        weights = don.weights
      )
    } else if (opts$type_predict == "loop") {
      # Do the prediction using a rcpp loop
      pred_raw <- statMatch.CCA.phase1.predict(
        CV_X_A = fit$CV_X_A, # canonical variable for data set A
        CV_X_B = fit$CV_X_B, # canonical variable for data set B
        h = tune$results$h[1],
        d = tune$results$d[1],
        Y = data$mtx.don.NCV.n,
        comp.mtx = comp.mtx,
        don.weights = don.weights,
        opts = opts
      )
      colnames(pred_raw) <- colnames(data$mtx.don.NCV.n)
    }
    pred <- mtx2df(pred_raw, var_names = names.NCV.cat)
    # Compute the sum of distances between the predictions and their compatible vectors in the donor data set
    #   This vector is used when averaging multiple models, for stability reasons.
    sum_dist_comp_donor <- compute_sum_dist_comp_donor(pred_raw, matrix(unlist(data$mtx.don.NCV.raw), nrow = nrow(comp.mtx)), comp.mtx)

    out <- list(
      tune = tune,
      fit = fit,
      df.match.don = pred,
      df.match = cbind(df.rec, pred),
      sum_dist_comp_donor = sum_dist_comp_donor
    )
  } else {
    if (length(names.NCV.cat) == 0) {
      message("\nPhase 1 - No categorical variables to predict. Passing to the next phase ...")
    }
    if (length(names.CV.cat) == 0 & length(names.NCV.cat) > 0) {
      stop("No common categorical variables. Impossible to compute compatibility matrix.")
    }
    out <- list(df.match = df.rec)
  }
  return(out)
}
