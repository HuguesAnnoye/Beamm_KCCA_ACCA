#' @title CCA - PHASE 2 TUNING AND PREDICTION in //
#' @description Tuning and prediction for CCA statistical matching process in //
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.CCA.options()} which contains all the options.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @return RETURN
#' @details DETAILS
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr %>%
#' @noRd

statMatch.CCA.phase2.tuning.par <- function(df.don, names.CV, names.NCV, don.weights, zero.constraints, opts) {

  # Prepare data for cross-validation
  input <- prepare_data_tuning(df.don, names.CV, names.NCV, don.weights, opts, method = "CCA")
  n_combs_new <- nrow(opts$P2$hpar)
  results <- opts$P2$hpar
  results$val_objfun <- NA
  results$time <- NA

  message("\nPhase 2. Tuning CCA ...")

  count_h <- 0
  fit <- list()
  # Perform cross-validation
  for (i in seq_len(n_combs_new)) {
    t_start <- Sys.time()
    # Count the number of h
    count_h <- count_h + 1
    # Do we fit the CCA
    if (count_h == 1) {
      fit_CCA <- TRUE
    } else {
      fit_CCA <- FALSE
    }

    # Predict and stack results over the k folds
    PRED.raw <- NULL
    cl <- makeCluster(opts$nc, outfile = "")
    # cl <- makeCluster(nc, setup_strategy = "sequential")#
    # clusterExport(cl, varlist = c("BEAMM.KCCAACCA.CCA.fit", "BEAMM.KCCAACCA.CCA.predict", "kccaw", "gauss3d", "gausskcca", "weighted.sd")) # , env=environment())
    registerDoParallel(cl)
    #for (fold in seq_len(opts$n_fold)) {
    res_fold <- foreach(fold = seq_len(opts$n_fold),.export=c("input","opts","fit", "fit_CCA"), .packages = c("BEAMM.KCCAACCA")) %dopar% {

      message(Sys.time(), ". Combination ", i, " of ", n_combs_new, ". Fold ", fold, " of ", opts$n_fold)
      if (fit_CCA) {
        fit_fold <- statMatch.CCA.fit(input$raw_data[[fold]]$mtx.don.CV.raw,
          input$raw_data[[fold]]$mtx.don.NCV.raw,
          weights = input$weights[[fold]],
          X2 = input$raw_data[[fold]]$mtx.rec.CV.raw,
          d = opts$P2$hpar$d[i]
        )
        fit_tuning <- fit_fold
      } else {
        # Compute predictions based on compatibility matrix
        fit_tuning <- fit[[fold]]
      }
      if (opts$type_predict == "matrix") {
        # Do the prediction using a matrix
        pred <- statMatch.KCCA.phase2.predict(
          CV_X_A = fit_tuning$CV_X_A, # canonical variables for data set A
          CV_X_B = fit_tuning$CV_X_B, # canonical variables for data set B
          Y = input$raw_data[[fold]]$mtx.don.NCV.n,
          scaling = opts$scaling,
          h = opts$P2$hpar$h[i],
          d = opts$P2$hpar$d[i],
          kernel_predict = opts$P2$kernel_predict_tuning,
          rot = opts$rot,
          matrix.tot.possibilitities = NULL, # input$comp.idx[[fold]],
          names.NCV = colnames(input$raw_data[[fold]]$mtx.don.NCV.n),
          weights = input$weights[[fold]],
          data_cat_rec = mtx2df(input$raw_data[[fold]]$mtx.rec.CV.n, names.CV),
          data_cat_don = mtx2df(input$raw_data[[fold]]$mtx.don.CV.n, names.CV),
          zero.constraints = zero.constraints,
          print.details = opts$print.details_tuning
        )
      } else if (opts$type_predict == "loop") {
        # Do the prediction using a rcpp loop
        pred <- statMatch.CCA.phase2.predict(
          CV_X_A = fit_tuning$CV_X_A, # canonical variables for data set A
          CV_X_B = fit_tuning$CV_X_B, # canonical variables for data set B
          h = opts$P2$hpar$h[i],
          d = opts$P2$hpar$d[i],
          mtx.rec.CV.raw = input$raw_data[[fold]]$mtx.rec.CV.n,
          mtx.don.CV.raw = input$raw_data[[fold]]$mtx.don.CV.n,
          mtx.don.NCV.raw = input$raw_data[[fold]]$mtx.don.NCV.n,
          don.weights = input$weights[[fold]],
          opts = opts,
          zero.constraints = zero.constraints,
          keep_unconstrained = F
        )
        colnames(pred) <- colnames(input$raw_data[[fold]]$mtx.don.NCV.n)
      }
      # Stack the predictions
      #PRED.raw <- rbind(PRED.raw, pred)
      return(list(pred=pred,fit=fit_tuning))
    }
    stopCluster(cl)
    for (fold in seq_len(opts$n_fold)) {
      PRED.raw <- rbind(PRED.raw,res_fold[[fold]]$pred)
      fit[[fold]] <- res_fold[[fold]]$fit
    }
    res_fold <- NULL
    # Store elapsed time
    results$time[i] <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    # Compute the objective function
    results$val_objfun[i] <- compute_obj_function(input$mtx.don.NCV.n, PRED.raw[, colnames(input$mtx.don.NCV.n)], don.weights, type = opts$P2$objmethod)
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
