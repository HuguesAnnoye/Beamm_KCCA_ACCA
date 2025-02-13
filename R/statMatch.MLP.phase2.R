#' @title statMatch.MLP.phase2
#' @description Statistical matching (phase2) using MMLP
#' @param df.rec a data.frame, the receiver data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param opts a list returned by the function \code{statMatch.MLP.options()} which contains all the options.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#'
#' @importFrom dplyr select all_of
#' @importFrom magrittr '%>%'
#' @importFrom keras k_clear_session fit
#' @importFrom tensorflow set_random_seed
#' @noRd

statMatch.MLP.phase2 <- function(df.rec, df.don, don.weights, names.CV, names.NCV, zero.constraints, opts) {

  if (opts$print.details) {
    message("\nPhase 2 - Matching continuous variables.")
    message("\nCommon variables: ", paste0(names.CV, collapse = ", "), ".")
    message("\nNon-common variables: ", paste0(names.NCV, collapse = ", "), ".")
  }

  # Tune the hyperparameters (if necessary)
  if (opts$P2$n_combs != 1)
    tune <- statMatch.MLP.phase2.tuning(df.don, names.CV, names.NCV, don.weights, zero.constraints, opts)
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
  set_random_seed(opts$seed_phase2)
  # Initialize the MMLP with keras
  model_phase2 <- statMatch.MLP.phase2.model(
    n_X = ncol(data$mtx.don.CV.raw),
    n_Y = ncol(data$mtx.don.NCV.raw),
    n_hidden_layers = opts$P2$nlayers,
    n_units = as.numeric((tune$results %>% select(starts_with("units")))[1,]),
    penL1 = tune$results$penL1[1],
    penL2 = tune$results$penL2[1],
    learning_rate = tune$results$lr[1]
  )
  # Fit model with keras
  fit <- model_phase2 %>%
    fit(
      x = data$mtx.don.CV.raw,
      y = data$mtx.don.NCV.raw,
      batch_size = opts$P2$batch_size,
      epochs = opts$P2$epochs,
      # validation_data = list(
      #   data$mtx.don.CV.raw,
      #   data$mtx.don.NCV.raw,
      #   matrix(don.weights)
      # ),
      verbose = 0,
      callbacks = opts$callbacks_keras,
      sample_weight = matrix(don.weights)
    )
  pred_raw <- statMatch.MLP.phase2.predict(model_phase2,
                                           data$mtx.rec.CV.raw,
                                           data$mtx.don.NCV.raw,
                                           don.weights,
                                           opts,
                                           zero.constraints,
                                           keep_unconstrained = T)
  pred <- mtx2df(pred_raw, var_names = names.NCV)

  out <- list(tune = tune,
              fit = fit,
              df.match.don = pred,
              df.match = cbind(df.rec, pred))
}
