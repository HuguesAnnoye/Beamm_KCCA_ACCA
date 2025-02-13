#' @title ACCA - PHASE 1 TUNING AND PREDICTION
#' @description Tuning and prediction for autoencoder and CCA statistical matching process
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.ACCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr %>%
#' @importFrom keras k_clear_session
#' @importFrom tensorflow set_random_seed
#' @noRd
#'
statMatch.ACCA.phase1.tuning <- function(df.don, names.CV, names.NCV, don.weights, opts) {
  # Prepare data for cross-validation
  names.CV.cat <- names.CV[unlist(lapply(df.don %>% select(all_of(names.CV)), is.factor))]
  input <- prepare_data_tuning(df.don, names.CV, names.NCV, don.weights, opts,
                               comp.mtx = T, names.CV.cat = names.CV.cat)
  n_X <- ncol(input$raw_data[[1]]$mtx.don.CV.raw)
  n_Y <- ncol(input$raw_data[[1]]$mtx.don.NCV.raw)

  n_combs_new <- nrow(opts$P1$hpar)
  #autoencoder_X <- vector(mode = "list", opts$n_fold)
  #autoencoder_Y <- vector(mode = "list", opts$n_fold)
  cca_model <- vector(mode = "list", opts$n_fold)
  results <- opts$P1$hpar
  results$val_objfun <- NA
  results$avg_epochs_X <- NA
  results$avg_epochs_Y <- NA
  results$time <- NA
  info <- list(
    epochs_X = matrix(NA, n_combs_new, opts$n_fold),
    epochs_Y = matrix(NA, n_combs_new, opts$n_fold)
  )

  message("\nPhase 1. Tuning ACCA ...")
  count_h <- 0
  # Perform cross-validation
  for (i in seq_len(n_combs_new)) {
    t_start <- Sys.time()
    # Count the number of h
    count_h <- count_h + 1
    if (count_h == opts$P1$n_h + 1) count_h <- 1
    # Do we fit the autoencoder ?
    if (count_h == 1) fit_autoencoder <- TRUE
    else fit_autoencoder <- FALSE

    # Predict and stack results over the k folds
    PRED.raw <- NULL
    for (fold in seq_len(opts$n_fold)) {
      message(Sys.time(), ". Combination ", i, " of ", n_combs_new, ". Fold ", fold, " of ", opts$n_fold)
      if (fit_autoencoder) {
        # Clear previous keras session
        k_clear_session()
        # Set the same seed for each combination of hyperparameters
        if (fold == 1) set_random_seed(opts$seed_phase1)
        # Initialize the autoencoder for the common variables
        autoencoder_X <- statMatch.ACCA.autoencoder(
          n_var = n_X,
          n_lat = opts$P1$hpar$lat_X[i],
          n_hidden_layers = opts$P1$nlayers,
          n_units = as.numeric((opts$P1$hpar %>% select(starts_with("units_X_")))[i,]),
          penL1 = opts$P1$hpar$penL1_X[i],
          penL2 = opts$P1$hpar$penL2_X[i],
          lr = opts$P1$hpar$lr_X[i]
        )
        # Initialize the autoencoder for the non-common variables
        autoencoder_Y <- statMatch.ACCA.autoencoder(
          n_var = n_Y,
          n_lat = opts$P1$hpar$lat_Y[i],
          n_hidden_layers = opts$P1$nlayers,
          n_units = as.numeric((opts$P1$hpar %>% select(starts_with("units_Y_")))[i,]),
          penL1 = opts$P1$hpar$penL1_Y[i],
          penL2 = opts$P1$hpar$penL2_Y[i],
          lr = opts$P1$hpar$lr_Y[i]
        )
        # Fit the autoencoder for the common variables
        history_X <- autoencoder_X$model %>%
          fit(
            x = input$raw_data[[fold]]$mtx.don.CV.raw,
            y = input$raw_data[[fold]]$mtx.don.CV.raw,
            batch_size = opts$P1$batch_size,
            epochs = opts$P1$epochs,
            # validation_data = list(
            #   input$raw_data[[fold]]$mtx.don.CV.raw,
            #   input$raw_data[[fold]]$mtx.don.CV.raw,
            #   matrix(input$weights[[fold]])
            # ),
            verbose = 0,
            callbacks = opts$callbacks_keras,
            sample_weight = matrix(input$weights[[fold]])
          )
        # Fit the autoencoder for the non-common variables
        history_Y <- autoencoder_Y$model %>%
          fit(
            x = input$raw_data[[fold]]$mtx.don.NCV.raw,
            y = input$raw_data[[fold]]$mtx.don.NCV.raw,
            batch_size = opts$P1$batch_size,
            epochs = opts$P1$epochs,
            # validation_data = list(
            #   input$raw_data[[fold]]$mtx.don.CV.raw,
            #   input$raw_data[[fold]]$mtx.don.CV.raw,
            #   matrix(input$weights[[fold]])
            # ),
            verbose = 0,
            callbacks = opts$callbacks_keras,
            sample_weight = matrix(input$weights[[fold]])
          )
        # Store the number of epochs
        info$epochs_X[i, fold] <- length(history_X$metrics$loss)
        info$epochs_Y[i, fold] <- length(history_Y$metrics$loss)

        cca_model[[fold]] <- statMatch.ACCA.phase1.cca(autoencoder_X,
                                               autoencoder_Y,
                                               input$raw_data[[fold]]$mtx.rec.CV.raw,
                                               input$raw_data[[fold]]$mtx.don.CV.raw,
                                               input$raw_data[[fold]]$mtx.don.NCV.raw,
                                               input$weights[[fold]],
                                               d=opts$P1$hpar$d[i],
                                               opts)

      } else { message(cca_model[[fold]]$can_var_CCA_mess) }
      # Compute predictions based on compatibility matrix

      pred <- statMatch.ACCA.phase1.predict(cca_model[[fold]],
                                            opts$P1$hpar$h[i],
                                            d=opts$P1$hpar$d[i],
                                            input$raw_data[[fold]]$mtx.don.NCV.raw,
                                            input$comp.idx[[fold]],
                                            input$weights[[fold]],
                                            opts)
      # Stack the predictions
      PRED.raw <- rbind(PRED.raw, pred)
    }
    # Store elapsed time
    results$time[i] <-  as.numeric(difftime(Sys.time(), t_start, units = 'mins'))
    # Compute average of the epochs over the k folds
    results$avg_epochs_X[i] <- mean(info$epochs_X[i,])
    results$avg_epochs_Y[i] <- mean(info$epochs_Y[i,])
    # Compute the objective function
    results$val_objfun[i] <- compute_obj_function(input$mtx.don.NCV.raw, PRED.raw, don.weights, type = "wRMSE")
    if (opts$print.details) print(results[i,])
  }

  # Order results and info$epochs by val_objfun in increasing order
  idx_objfun <- order(results$val_objfun)
  results <- results[idx_objfun,]
  info$epochs_X <- info$epochs_X[idx_objfun,]
  info$epochs_Y <- info$epochs_Y[idx_objfun,]

  out <- list(
    results = results,
    info = info
  )

  return(out)
}
