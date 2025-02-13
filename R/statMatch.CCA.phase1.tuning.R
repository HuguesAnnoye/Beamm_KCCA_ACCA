#' @title CCA - PHASE 1 TUNING
#' @description Tuning and prediction for CCA statistical matching process
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.CCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr %>%
#' @importFrom keras k_clear_session
#' @importFrom tensorflow set_random_seed
#' @noRd

statMatch.CCA.phase1.tuning <- function(df.don, names.CV, names.NCV, don.weights, opts) {

  # Prepare data for cross-validation
  names.CV.cat <- names.CV[unlist(lapply(df.don %>% select(all_of(names.CV)), is.factor))]
  input <- prepare_data_tuning(df.don, names.CV, names.NCV, don.weights, opts,
    comp.mtx = T, names.CV.cat = names.CV.cat, method = "CCA"
  )

  n_combs_new <- nrow(opts$P1$hpar)
  results <- opts$P1$hpar
  results$val_objfun <- NA
  results$time <- NA

  message("\nPhase 1. Tuning CCA ...")

  count_h <- 0
  # Perform cross-validation
  fit <- list()
  for (i in seq_len(n_combs_new)) {
    t_start <- Sys.time()
    # Count the number of h
    count_h <- count_h + 1
    # Do we fit the autoencoder ?
    if (count_h == 1) {
      fit_CCA <- TRUE
    } else {
      fit_CCA <- FALSE
    }

    # Predict and stack results over the k folds
    PRED.raw <- NULL
    for (fold in seq_len(opts$n_fold)) {
      message(Sys.time(), ". Combination ", i, " of ", n_combs_new, ". Fold ", fold, " of ", opts$n_fold)
      if (fit_CCA) {
        fit[[fold]] <- statMatch.CCA.fit(input$raw_data[[fold]]$mtx.don.CV.raw,
          input$raw_data[[fold]]$mtx.don.NCV.raw,
          weights = input$weights[[fold]],
          X2 = input$raw_data[[fold]]$mtx.rec.CV.raw,
          d = opts$P1$hpar$d[i]
        )
      }
      # Compute predictions based on compatibility matrix
      fit_tuning <- fit[[fold]]
      if (opts$type_predict == "matrix") {
        # Do the prediction using a matrix
        pred <- statMatch.KCCA.phase1.predict(
          CV_X_A = fit_tuning$CV_X_A, # canonical variables for data set A
          CV_X_B = fit_tuning$CV_X_B, # canonical variables for data set B
          Y = input$raw_data[[fold]]$mtx.don.NCV.n,
          scaling = opts$scaling,
          h = opts$P1$hpar$h[i],
          d = opts$P1$hpar$d[i],
          kernel_predict = opts$P1$kernel_predict_tuning,
          rot = opts$rot,
          matrix.tot.possibilitities = input$comp.idx[[fold]],
          names.NCV = colnames(input$raw_data[[fold]]$mtx.don.NCV.n),
          weights = input$weights[[fold]]
        )
      } else if (opts$type_predict == "loop") {
        # Do the prediction using a rcpp loop
        pred <- statMatch.CCA.phase1.predict(
          CV_X_A = fit_tuning$CV_X_A, # canonical variable for data set A
          CV_X_B = fit_tuning$CV_X_B, # canonical variable for data set B
          h = opts$P1$hpar$h[i],
          d = opts$P1$hpar$d[i],
          Y = input$raw_data[[fold]]$mtx.don.NCV.n,
          comp.mtx = input$comp.idx[[fold]],
          don.weights = input$weights[[fold]],
          opts = opts
        )
        colnames(pred) <- colnames(input$raw_data[[fold]]$mtx.don.NCV.n)
      }

      # Stack the predictions
      PRED.raw <- rbind(PRED.raw, pred)
    }
    # Store elapsed time
    results$time[i] <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    # Compute the objective function
    if (opts$P1$objmethod == "wMCR") {
      val0 <- rep(NA, length(colnames(input$mtx.don.NCV.n)))
      for (i2 in 1:length(colnames(input$mtx.don.NCV.n))) {
        nomi <- colnames(input$mtx.don.NCV.n)[i2]
        val0[i2] <- compute_MCR(input$mtx.don.NCV.n[, nomi], PRED.raw[, nomi], weights = don.weights)
      }
      results$val_objfun[i] <- mean(val0)
    } else {
      results$val_objfun[i] <- compute_obj_function(input$mtx.don.NCV.n, PRED.raw[, colnames(input$mtx.don.NCV.n)], don.weights, type = opts$P1$objmethod)
    }
    if (opts$print.details) print(results[i, ])
  }

  # Order results and info$epochs by val_objfun in increasing order
  idx_objfun <- order(results$val_objfun)
  results <- results[idx_objfun, ]

  out <- list(
    results = results,
    idx_objfun = idx_objfun
  )

  return(out)
}
