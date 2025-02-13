#' @title PCA - PHASE 2 TUNING AND PREDICTION
#' @description Tuning and prediction for PCA - CCA statistical matching process
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param don.weights a numeric vector, the sample weights of tkhe donor data set.
#' @param opts a list returned by the function \code{statMatch.PCA.options()} which contains all the options.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @return RETURN
#' @details DETAILS
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr %>%
#' @noRd

statMatch.PCA.phase2.tuning <- function(df.don, names.CV, names.NCV, don.weights, zero.constraints, opts) {

  # Prepare data for cross-validation
  input <- prepare_data_tuning(df.don, names.CV, names.NCV, don.weights, opts)

  n_combs_new <- nrow(opts$P2$hpar)
  results <- opts$P2$hpar
  results$val_objfun <- NA
  results$time <- NA

  message("\nPhase 2. Tuning PCA ...")

  count_h <- 0
  fit <- list()
  # Perform cross-validation
  for (i in seq_len(n_combs_new)) {
    t_start <- Sys.time()
    # Count the number of h
    count_h <- count_h + 1
    if (count_h == opts$P2$n_h + 1) count_h <- 1
    # Do we fit the PCA ?
    if (count_h == 1) {
      fit_PCA <- TRUE
    } else {
      fit_PCA <- FALSE
    }

    # Predict and stack results over the k folds
    PRED.raw <- NULL
    for (fold in seq_len(opts$n_fold)) {
      message(Sys.time(), ". Combination ", i, " of ", n_combs_new, ". Fold ", fold, " of ", opts$n_fold)
      if (fit_PCA) {
        fit[[fold]] <- statMatch.PCA.fit(input$raw_data[[fold]]$mtx.don.CV.raw,
          input$raw_data[[fold]]$mtx.don.NCV.raw,
          weights = input$weights[[fold]],
          X2 = input$raw_data[[fold]]$mtx.rec.CV.raw,
          d = opts$P2$hpar$d[i],
          lat_X = opts$P2$hpar$lat_X[i],
          lat_Y = opts$P2$hpar$lat_Y[i]
        )
      }
      # Compute predictions based on compatibility matrix
      fit_tuning <- fit[[fold]]

      if (opts$type_predict == "matrix") {
        # Do the prediction using a matrix
        # Go back to the original scale of the variables used to create the compatibility matrix for ZC variables
        if (opts$scaling == "z-score") {
          data_cat_rec <- scale_mtx_inv(input$raw_data[[fold]]$mtx.rec.CV.raw,
            center = attr(input$raw_data[[fold]]$mtx.rec.CV.raw, "scaled:center"),
            scale = attr(input$raw_data[[fold]]$mtx.rec.CV.raw, "scaled:scale")
          )
          data_cat_don <- scale_mtx_inv(input$raw_data[[fold]]$mtx.don.CV.raw,
            center = attr(input$raw_data[[fold]]$mtx.don.CV.raw, "scaled:center"),
            scale = attr(input$raw_data[[fold]]$mtx.don.CV.raw, "scaled:scale")
          )
        } else if (opts$scaling == "min-max") {
          data_cat_rec <- normalize_mtx_inv(input$raw_data[[fold]]$mtx.rec.CV.raw,
            min = attr(input$raw_data[[fold]]$mtx.rec.CV.raw, "min"),
            max = attr(input$raw_data[[fold]]$mtx.rec.CV.raw, "max")
          )
          data_cat_don <- normalize_mtx_inv(input$raw_data[[fold]]$mtx.don.CV.raw,
            min = attr(input$raw_data[[fold]]$mtx.don.CV.raw, "min"),
            max = attr(input$raw_data[[fold]]$mtx.don.CV.raw, "max")
          )
        } else {
          data_cat_rec <- input$raw_data[[fold]]$mtx.rec.CV.raw
          data_cat_don <- input$raw_data[[fold]]$mtx.don.CV.raw
        }
        data_cat_rec <- mtx2df(data_cat_rec, names.CV)
        data_cat_don <- mtx2df(data_cat_don, names.CV)
        pred <- statMatch.KCCA.phase2.predict(
          CV_X_A = fit_tuning$CV_X_A, # canonical variables for data set A
          CV_X_B = fit_tuning$CV_X_B, # canonical variables for data set B
          Y = input$raw_data[[fold]]$mtx.don.NCV.raw,
          h = opts$P2$hpar$h[i],
          d = opts$P2$hpar$d[i],
          kernel_predict = opts$P2$kernel_predict_tuning,
          scaling = opts$scaling,
          rot = opts$rot,
          matrix.tot.possibilitities = NULL, # input$comp.idx[[fold]],
          names.NCV = colnames(input$raw_data[[fold]]$mtx.don.NCV.raw),
          weights = input$weights[[fold]],
          data_cat_rec = data_cat_rec,
          data_cat_don = data_cat_don,
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
          mtx.rec.CV.raw = input$raw_data[[fold]]$mtx.rec.CV.raw,
          mtx.don.CV.raw = input$raw_data[[fold]]$mtx.don.CV.raw,
          mtx.don.NCV.raw = input$raw_data[[fold]]$mtx.don.NCV.raw,
          don.weights = input$weights[[fold]],
          opts = opts,
          zero.constraints = zero.constraints,
          keep_unconstrained = F
        )
        colnames(pred) <- colnames(input$raw_data[[fold]]$mtx.don.NCV.raw)
      }
      # Stack the predictions
      PRED.raw <- rbind(PRED.raw, pred)
    }
    # Store elapsed time
    results$time[i] <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    # Compute the objective function
    results$val_objfun[i] <- compute_obj_function(input$mtx.don.NCV.raw, PRED.raw[, colnames(input$mtx.don.NCV.raw)], don.weights, type = opts$P2$objmethod)
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
