#' @title statMatch.MLP.phase1.tuning
#' @description Tuning (phase1) of the statistical matching algorithm using MLP
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param opts a list, supplied using the function \code{statMatch.MLP.options()}, containing all the options.
#' @return A data.frame, ...
#'
#' @importFrom magrittr %>%
#' @importFrom keras k_clear_session fit
#' @importFrom tensorflow set_random_seed
#' @noRd

statMatch.MLP.phase1.tuning <- function(df.don, names.CV, names.NCV, don.weights, opts) {

  # Prepare data for cross-validation
  names.CV.cat <- names.CV[unlist(lapply(df.don %>% select(all_of(names.CV)), is.factor))]
  input <- prepare_data_tuning(df.don, names.CV, names.NCV, don.weights, opts, method = "MLP_phase1",
                               comp.mtx = T, names.CV.cat = names.CV.cat)
  n_X <- ncol(input$raw_data[[1]]$mtx.don.CV.raw)
  n_Y <- length(names.NCV)
  n_levels_Y <- sapply(input$raw_data[[1]]$mtx.don.NCV.raw, ncol)

  results <- opts$P1$hpar
  results$val_objfun <- NA
  results$avg_epochs <- NA
  results$time <- NA
  info <- list(epochs = matrix(NA, opts$P1$n_combs, opts$n_fold))

  message("\nPhase 1. Tuning MLP ...")

  # Perform cross-validation
  for (i in seq_len(opts$P1$n_combs)) {
    t_start <- Sys.time()
    # Clear previous keras session
    k_clear_session()
    # Set the same seed for each combination of hyperparameters
    set_random_seed(opts$seed_phase1)
    # Initialize the MMLP with keras
    model <- statMatch.MLP.phase1.model(
      n_X, n_Y, n_levels_Y,
      n_hidden_layers = opts$P1$nlayers,
      n_units = as.numeric((opts$P1$hpar %>% select(starts_with("units")))[i,]),
      penL1 = opts$P1$hpar$penL1[i],
      penL2 = opts$P1$hpar$penL2[i],
      learning_rate = opts$P1$hpar$lr[i]
    )

    # Predict and stack results over the k folds
    PRED.raw <- NULL
    for (fold in seq_len(opts$n_fold)) {
      message(Sys.time(), ". Combination ", i, " of ", opts$P1$n_combs, ". Fold ", fold, " of ", opts$n_fold)
      # Fit model with keras
      history <- model %>%
        fit(
          x = input$raw_data[[fold]]$mtx.don.CV.raw,
          y = input$raw_data[[fold]]$mtx.don.NCV.raw,
          batch_size = opts$P1$batch_size,
          epochs = opts$P1$epochs,
          # validation_data = list(
          #   input$raw_data[[fold]]$mtx.don.CV.raw,
          #   input$raw_data[[fold]]$mtx.don.NCV.raw,
          #   matrix(input$weights[[fold]])
          # ),
          verbose = 0,
          callbacks = opts$callbacks_keras,
          sample_weight = matrix(input$weights[[fold]])
        )
      info$epochs[i, fold] <- length(history$metrics$loss)
      # Compute predictions based on compatibility matrix
      pred <- statMatch.MLP.phase1.predict(model,
                                           input$raw_data[[fold]]$mtx.rec.CV.raw,
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
    results$avg_epochs[i] <- mean(info$epochs[i,])
    # Compute the objective function
    results$val_objfun[i] <- compute_obj_function(input$mtx.don.NCV.raw, PRED.raw, don.weights, type = "wRMSE")
    if (opts$print.details) print(results[i,])
  }

  # Order results and info$epochs by val_objfun in increasing order
  idx_objfun <- order(results$val_objfun)
  results <- results[idx_objfun,]
  info$epochs <- info$epochs[idx_objfun,]

  out <- list(
    results = results,
    info = info
  )

  return(out)
}
