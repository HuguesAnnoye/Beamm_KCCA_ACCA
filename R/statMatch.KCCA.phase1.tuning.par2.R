#' @title KCCA - PHASE 1 TUNING AND PREDICTION
#' @description Tuning and prediction for Kernel CCA statistical matching process
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.KCCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr %>%
#' @noRd

statMatch.KCCA.phase1.tuning.par2 <- function(df.don, names.CV, names.NCV, don.weights, opts) {

  # Prepare data for cross-validation
  names.CV.cat <- names.CV[unlist(lapply(df.don %>% select(all_of(names.CV)), is.factor))]
  input <- prepare_data_tuning(
    df.don, names.CV, names.NCV, don.weights, opts,
    comp.mtx = T, names.CV.cat = names.CV.cat
  )
  n_combs_new <- nrow(opts$P1$hpar)
  results <- opts$P1$hpar
  results$val_objfun <- NA
  results$time <- NA

  message("\nPhase 1. Tuning KCCA ...")

  count_h <- 0
  # Perform cross-validation
  cl <- makeCluster(opts$nc, outfile = "")
  #cl <- makeCluster(opts$nc, setup_strategy = "sequential")
  registerDoParallel(cl)
  results_for <- foreach(
    j = seq_len(opts$P1$n_combs), #.export = c("count_h", "input", "opts", "results"),
    # c("input","opts","fit", "fit_KCCA"),
    .combine = rbind,
    .inorder = TRUE,
    .packages = c("BEAMM.KCCAACCA")
  ) %dopar% {
    fit <- list()
    results2 <- results[(opts$P1$n_h * (j - 1) + 1):(opts$P1$n_h * (j)), ]
    for (ih in seq_len(opts$P1$n_h)) {
      i <- ih + opts$P1$n_h * (j - 1)
      t_start <- Sys.time()
      # Count the number of h
      count_h <- count_h + 1
      if (count_h == opts$P1$n_h + 1) count_h <- 1
      # Do we fit the KCCA ?
      if (count_h == 1) {
        fit_KCCA <- TRUE
      } else {
        fit_KCCA <- FALSE
      }

      # Predict and stack results over the k folds
      PRED.raw <- NULL
      for (fold in seq_len(opts$n_fold)) {
        message(Sys.time(), ". Combination ", i, " of ", n_combs_new, ". Fold ", fold, " of ", opts$n_fold)
        if (fit_KCCA) {
          fit[[fold]] <- statMatch.KCCA.fit(input$raw_data[[fold]]$mtx.don.CV.raw,
            input$raw_data[[fold]]$mtx.don.NCV.raw,
            weights = input$weights[[fold]],
            X2 = input$raw_data[[fold]]$mtx.rec.CV.raw,
            d = opts$P1$hpar$d[i],
            h = opts$P1$hpar$hx[i],
            hy = opts$P1$hpar$hy[i],
            g = opts$P1$hpar$g[i],
            rot = opts$rot
          )
        }
        # Compute predictions based on compatibility matrix
        fit_tuning <- fit[[fold]]
        if (opts$type_predict == "matrix") {
          # Do the prediction using a matrix
          pred <- statMatch.KCCA.phase1.predict(
            CV_X_A = fit_tuning$CV_X_A, # canonical variables for data set A
            CV_X_B = fit_tuning$CV_X_B, # canonical variables for data set B
            Y = input$raw_data[[fold]]$mtx.don.NCV.raw,
            h = opts$P1$hpar$h[i],
            d = opts$P1$hpar$d[i],
            kernel_predict = opts$P1$kernel_predict_tuning,
            scaling = opts$scaling,
            rot = opts$rot,
            matrix.tot.possibilitities = input$comp.idx[[fold]],
            names.NCV = colnames(input$raw_data[[fold]]$mtx.don.NCV.raw),
            weights = input$weights[[fold]]
          )
        } else if (opts$type_predict == "loop") {
          # Do the prediction using a rcpp loop
          pred <- statMatch.CCA.phase1.predict(
            CV_X_A = fit_tuning$CV_X_A, # canonical variable for data set A
            CV_X_B = fit_tuning$CV_X_B, # canonical variable for data set B
            h = opts$P1$hpar$h[i],
            d = opts$P1$hpar$d[i],
            Y = input$raw_data[[fold]]$mtx.don.NCV.raw,
            comp.mtx = input$comp.idx[[fold]],
            don.weights = input$weights[[fold]],
            opts = opts
          )
          colnames(pred) <- colnames(input$raw_data[[fold]]$mtx.don.NCV.raw)
        }
        # Stack the predictions
        PRED.raw <- rbind(PRED.raw, pred)
        pred <- NULL
      }
      # Store elapsed time
      results2$time[ih] <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      # Compute the objective function
      if (opts$P1$objmethod == "wMCR") {
        val0 <- rep(NA, length(colnames(input$mtx.don.NCV.raw)))
        for (i2 in 1:length(colnames(input$mtx.don.NCV.raw))) {
          nomi <- colnames(input$mtx.don.NCV.raw)[i2]
          val0[i2] <- compute_MCR(input$mtx.don.NCV.raw[, nomi], PRED.raw[, nomi], weights = don.weights)
        }
        results2$val_objfun[ih] <- mean(val0)
      } else {
        results2$val_objfun[ih] <- compute_obj_function(input$mtx.don.NCV.raw, PRED.raw[, colnames(input$mtx.don.NCV.raw)], don.weights, type = opts$P1$objmethod)
      }
      if (opts$print.details) print(results2[ih, ])
      PRED.raw <- NULL
      }
    return(results2)
  }
  stopCluster(cl)
  # Order results and info$epochs by val_objfun in increasing order
  idx_objfun <- order(results_for$val_objfun)
  results <- results_for[idx_objfun, ]
  if (opts$print.details) {
    print(paste0("Best parameters for phase 1 is "))
    print(results[1, ])
  }
  out <- list(
    results = results,
    idx_objfun = idx_objfun
  )

  return(out)
}
