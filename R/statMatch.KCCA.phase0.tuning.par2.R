#' @title KCCA - phase 0 TUNING AND PREDICTION
#' @description Tuning and prediction for Kernel CCA statistical matching process
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param names.NCV_cont a character vector with the names of the non-common variables used for statistical matching.
#' @param names.NCV_cat a character vector with the names of the non-common variables used for statistical matching.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.KCCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr %>%
#' @noRd

statMatch.KCCA.phase0.tuning.par2 <- function(df.don, names.CV, names.NCV, names.NCV_cat, names.NCV_cont, don.weights, zero.constraints, opts) {

  # Prepare data for cross-validation
  names.CV.cat <- names.CV[unlist(lapply(df.don %>% select(all_of(names.CV)), is.factor))]
  input <- prepare_data_tuning_phase0(
    df.don, names.CV, names.NCV,
    names.NCV_cat, names.NCV_cont, don.weights, opts,
    comp.mtx = T, names.CV.cat = names.CV.cat
  )
  n_combs_new <- nrow(opts$P0$hpar)
  results <- opts$P0$hpar
  results$val_objfun <- NA
  results$time <- NA

  message("\nPhase 0. Tuning KCCA ...")

  count_h <- 0
  # Perform cross-validation
  cl <- makeCluster(opts$nc, outfile = "")
  # cl <- makeCluster(opts$nc, setup_strategy = "sequential")
  registerDoParallel(cl)
  results_for <- foreach(
    j = seq_len(opts$P0$n_combs), # .export = c("count_h", "input", "opts", "results"),
    # c("input","opts","fit", "fit_KCCA"),
    .combine = rbind,
    .inorder = TRUE,
    .packages = c("BEAMM.KCCAACCA")
  ) %dopar% {
    fit <- list()
    results2 <- results[(opts$P0$n_h * (j - 1) + 1):(opts$P0$n_h * (j)), ]
    for (ih in seq_len(opts$P0$n_h)) {
      i <- ih + opts$P0$n_h * (j - 1)
      t_start <- Sys.time()
      # Count the number of h
      count_h <- count_h + 1
      if (count_h == opts$P0$n_h + 1) count_h <- 1
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
            d = opts$P0$hpar$d[i],
            h = opts$P0$hpar$hx[i],
            hy = opts$P0$hpar$hy[i],
            g = opts$P0$hpar$g[i],
            rot = opts$rot
          )
        }
        # Compute predictions based on compatibility matrix
        fit_tuning <- fit[[fold]]
        if (opts$type_predict == "matrix") {
          # Do the prediction using a matrix
          pred_cat <- statMatch.KCCA.phase1.predict(
            CV_X_A = fit_tuning$CV_X_A, # canonical variables for data set A
            CV_X_B = fit_tuning$CV_X_B, # canonical variables for data set B
            Y = input$raw_data[[fold]]$mtx.don.NCV.raw_cat,
            h = opts$P0$hpar$h_cat[i],
            d = opts$P0$hpar$d[i],
            kernel_predict = opts$P0$kernel_predict_cat,
            scaling = opts$scaling,
            rot = opts$rot,
            matrix.tot.possibilitities = input$comp.idx[[fold]],
            names.NCV = colnames(input$raw_data[[fold]]$mtx.don.NCV.raw_cat),
            weights = input$weights[[fold]]
          )
          # pred_cat<- mtx2df(pred_cat, var_names = names.NCV_cat)
        } else if (opts$type_predict == "loop") {
          # Do the prediction using a rcpp loop
          pred_cat <- statMatch.CCA.phase1.predict(
            CV_X_A = fit_tuning$CV_X_A, # canonical variable for data set A
            CV_X_B = fit_tuning$CV_X_B, # canonical variable for data set B
            h = opts$P0$hpar$h_cat[i],
            Y = input$raw_data[[fold]]$mtx.don.NCV.raw_cat,
            comp.mtx = input$comp.idx[[fold]],
            don.weights = input$weights[[fold]],
            opts = opts
          )
          colnames(pred_cat) <- colnames(input$raw_data[[fold]]$mtx.don.NCV.raw_cat)
          # pred_cat<- mtx2df(pred_raw, var_names = names.NCV_cat)
        }

        # à modifier ; +Il faut que ce soit correct :

        input$raw_data[[fold]]$mtx.rec.CV.raw_cont <- cbind(input$raw_data[[fold]]$mtx.rec.CV.raw2, pred_cat) # df2mtx(cbind(data_cat_rec, pred))
        if (opts$scaling == "z-score") {
          input$raw_data[[fold]]$mtx.rec.CV.raw_cont <- scale_mtx(input$raw_data[[fold]]$mtx.rec.CV.raw_cont,
            center = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "scaled:center"),
            scale = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "scaled:scale")
          )
        } else if (opts$scaling == "min-max") {
          input$raw_data[[fold]]$mtx.rec.CV.raw_cont <- normalize_mtx(input$raw_data[[fold]]$mtx.rec.CV.raw_cont,
            min = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "min"),
            max = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "max")
          )
        }
        if (opts$type_predict == "matrix") {
          if (opts$scaling == "z-score") {
            data_cat_rec <- scale_mtx_inv(input$raw_data[[fold]]$mtx.rec.CV.raw_cont,
              center = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "scaled:center"),
              scale = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "scaled:scale")
            )
            data_cat_don <- scale_mtx_inv(input$raw_data[[fold]]$mtx.don.CV.raw_cont,
              center = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "scaled:center"),
              scale = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "scaled:scale")
            )
          } else if (opts$scaling == "min-max") {
            data_cat_rec <- normalize_mtx_inv(input$raw_data[[fold]]$mtx.rec.CV.raw_cont,
              min = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "min"),
              max = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "max")
            )
            data_cat_don <- normalize_mtx_inv(input$raw_data[[fold]]$mtx.don.CV.raw_cont,
              min = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "min"),
              max = attr(input$raw_data[[fold]]$mtx.don.CV.raw_cont, "max")
            )
          } else {
            data_cat_rec <- input$raw_data[[fold]]$mtx.rec.CV.raw_cont
            data_cat_don <- input$raw_data[[fold]]$mtx.don.CV.raw_cont
          }
          data_cat_rec <- mtx2df(data_cat_rec, c(names.CV, names.NCV_cat))
          data_cat_don <- mtx2df(data_cat_don, c(names.CV, names.NCV_cat))
          message(" Phase 2 for Combination ", i, " of ", n_combs_new, ". Fold ", fold, " of ", opts$n_fold)

          pred_cont <- statMatch.KCCA.phase2.predict(
            CV_X_A = fit_tuning$CV_X_A, # canonical variables for data set A
            CV_X_B = fit_tuning$CV_X_B, # canonical variables for data set B
            Y = input$raw_data[[fold]]$mtx.don.NCV.raw_cont,
            h = opts$P0$hpar$h[i],
            d = opts$P0$hpar$d[i],
            kernel_predict = opts$P0$kernel_predict,
            scaling = opts$scaling,
            rot = opts$rot,
            matrix.tot.possibilitities = NULL, # input$comp.idx[[fold]],
            names.NCV = colnames(input$raw_data[[fold]]$mtx.don.NCV.raw_cont),
            weights = input$weights[[fold]],
            data_cat_rec = data_cat_rec,
            data_cat_don = data_cat_don,
            zero.constraints = zero.constraints,
            print.details = opts$print.details_tuning
          )
        } else if (opts$type_predict == "loop") {
          pred_cont <- statMatch.CCA.phase2.predict(
            CV_X_A = fit_tuning$CV_X_A, # canonical variables for data set A
            CV_X_B = fit_tuning$CV_X_B, # canonical variables for data set B
            h = opts$P0$hpar$h[i],
            mtx.rec.CV.raw = input$raw_data[[fold]]$mtx.rec.CV.raw_cont,
            mtx.don.CV.raw = input$raw_data[[fold]]$mtx.don.CV.raw_cont,
            mtx.don.NCV.raw = input$raw_data[[fold]]$mtx.don.NCV.raw_cont,
            don.weights = input$weights[[fold]],
            opts = opts,
            zero.constraints = zero.constraints,
            keep_unconstrained = F
          )
          colnames(pred_cont) <- colnames(input$raw_data[[fold]]$mtx.don.NCV.raw_cont)
        }
        # pred_cont<- mtx2df(pred_raw, var_names = names.NCV_cont)
        # Stack the predictions
        PRED.raw <- rbind(PRED.raw, cbind(pred_cat, pred_cont))
        # voir si pred_cont doit pas être modifier
        pred_cat <- NULL
        pred_cont <- NULL
      }
      # Store elapsed time
      results2$time[ih] <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      # Compute the objective function
      if (opts$P0$objmethod == "wMCR") {
        val0 <- rep(NA, length(colnames(input$mtx.don.NCV.raw)))
        for (i2 in 1:length(colnames(input$mtx.don.NCV.raw))) {
          nomi <- colnames(input$mtx.don.NCV.raw)[i2]
          val0[i2] <- compute_MCR(input$mtx.don.NCV.raw[, nomi], PRED.raw[, nomi], weights = don.weights)
        }
        results2$val_objfun[ih] <- mean(val0)
      } else {
        message(opts$P0$objmethod, " for Combination ", i, " of ", n_combs_new, ". Fold ", fold, " of ", opts$n_fold)

        results2$val_objfun[ih] <- compute_obj_function(input$mtx.don.NCV.raw, PRED.raw[, colnames(input$mtx.don.NCV.raw)], don.weights, type = opts$P0$objmethod)
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
    print(paste0("Best parameters for phase 0 is "))
    print(results[1, ])
  }
  out <- list(
    results = results,
    idx_objfun = idx_objfun
  )

  return(out)
}
