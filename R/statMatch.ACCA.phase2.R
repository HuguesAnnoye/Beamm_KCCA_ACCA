#' @title statMatch.ACCA.phase2
#' @description Statistical matching (phase2) using ACCA
#' @param df.rec a data.frame, the receiver data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param opts a list returned by the function \code{statMatch.ACCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#'
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr '%>%'
#' @importFrom keras k_clear_session fit
#' @importFrom tensorflow set_random_seed
#' @noRd
#'
statMatch.ACCA.phase2 <- function(df.rec, df.don, don.weights, names.CV, names.NCV, zero.constraints, opts) {

  if (opts$print.details) {
    message("\nPhase 2 - Matching continuous variables.")
    message("\nCommon variables: ", paste0(names.CV, collapse = ", "), ".")
    message("\nNon-common variables: ", paste0(names.NCV, collapse = ", "), ".")
  }

  # Tune the hyperparameters (if necessary)
  if (opts$P2$n_combs != 1)
    tune <- statMatch.ACCA.phase2.tuning(df.don, names.CV, names.NCV, don.weights, zero.constraints, opts)
  else tune <- list(results = opts$P2$hpar)

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
  # Clear previous keras session
  k_clear_session()
  # Set same seed used in cross-validation
  set_random_seed(opts$seed_phase1)
  # Initialize the autoencoders with keras
  autoencoder_X <- statMatch.ACCA.autoencoder(
    n_var = ncol(data$mtx.don.CV.raw),
    n_lat = tune$results$lat_X[1],
    n_hidden_layers = opts$P2$nlayers,
    n_units = as.numeric((tune$results %>% select(starts_with("units_X_")))[1,]),
    penL1 = tune$results$penL1_X[1],
    penL2 = tune$results$penL2_X[1],
    lr = tune$results$lr_X[1]
  )
  autoencoder_Y <- statMatch.ACCA.autoencoder(
    n_var = ncol(data$mtx.don.NCV.raw),
    n_lat = tune$results$lat_Y[1],
    n_hidden_layers = opts$P2$nlayers,
    n_units = as.numeric((tune$results %>% select(starts_with("units_Y_")))[1,]),
    penL1 = tune$results$penL1_Y[1],
    penL2 = tune$results$penL2_Y[1],
    lr = tune$results$lr_Y[1]
  )
  fit_X <- autoencoder_X$model %>%
    fit(
      x = data$mtx.don.CV.raw,
      y = data$mtx.don.CV.raw,
      batch_size = opts$P2$batch_size,
      epochs = opts$P2$epochs,
      # validation_data = list(
      #   input$raw_data[[fold]]$mtx.don.CV.raw,
      #   input$raw_data[[fold]]$mtx.don.CV.raw,
      #   matrix(input$weights[[fold]])
      # ),
      verbose = 0,
      callbacks = opts$callbacks_keras,
      sample_weight = matrix(don.weights)
    )
  fit_Y <- autoencoder_Y$model %>%
    fit(
      x = data$mtx.don.NCV.raw,
      y = data$mtx.don.NCV.raw,
      batch_size = opts$P2$batch_size,
      epochs = opts$P2$epochs,
      # validation_data = list(
      #   input$raw_data[[fold]]$mtx.don.CV.raw,
      #   input$raw_data[[fold]]$mtx.don.CV.raw,
      #   matrix(input$weights[[fold]])
      # ),
      verbose = 0,
      callbacks = opts$callbacks_keras,
      sample_weight = matrix(don.weights)
    )
  # Compute CCA
  cca_raw <- statMatch.ACCA.phase2.cca(autoencoder_X, autoencoder_Y,
                                            data$mtx.rec.CV.raw,
                                            data$mtx.don.CV.raw,
                                            data$mtx.don.NCV.raw,
                                            don.weights,
                                            d=tune$results$d[1],
                                            opts, zero.constraints,
                                            keep_unconstrained = T)

  # Compute the predictions
  pred_raw <- statMatch.ACCA.phase2.predict(cca_raw,
                                            tune$results$h[1],
                                            tune$results$d[1],
                                            don.weights,
                                            opts, zero.constraints,
                                            keep_unconstrained = T)
  colnames(pred_raw) <- colnames(data$mtx.don.NCV.raw)
  pred <- mtx2df(pred_raw, var_names = names.NCV)
  out <- list(tune = tune,
              fit = fit,
              df.match.don = pred,
              df.match = cbind(df.rec, pred))
  return(out)
}
