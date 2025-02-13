#' @title ACCA - PHASE 1 TUNING AND PREDICTION
#' @description Tuning and prediction for autoencoder and CCA statistical matching process
#' @param df.rec a data.frame, the receiver data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param opts a list returned by the function \code{statMatch.ACCA.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr '%>%'
#' @importFrom keras k_clear_session
#' @importFrom tensorflow set_random_seed
#' @noRd
#'
statMatch.ACCA.phase1 <- function(df.rec, df.don, don.weights, names.CV, names.NCV, opts) {

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
    if (opts$P1$n_combs != 1)
      tune <- statMatch.ACCA.phase1.tuning(df.don, names.CV, names.NCV.cat, don.weights, opts)
    else tune <- list(results = opts$P1$hpar)

    if (opts$print.details)
      message("\nPhase 1. Fit and predict tuned model ...")

    # Prepare data
    data <- list()
    data$mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)))
    data$mtx.don.NCV.raw <- df2mtx(df.don %>% select(all_of(names.NCV.cat)))
    data$mtx.rec.CV.raw <- df2mtx(df.rec %>% select(all_of(names.CV)))

    if (opts$scaling != "no")
      data <- statMatch.scale(data, don.weights, opts$scaling)

    # Compute compatibility matrix
    comp.mtx <- compute_mtx_compatibilities(df.rec, df.don, names.CV.cat, opts$print.details, index = F)

    # Fit and predict using the best hyperparameters
    # Clear previous keras session
    k_clear_session()
    # Set same seed used in cross-validation
    set_random_seed(opts$seed_phase1)
    # Initialize the autoencoders with keras
    autoencoder_X <- statMatch.ACCA.autoencoder(
      n_var = ncol(data$mtx.don.CV.raw),
      n_lat = tune$results$lat_X[1],
      n_hidden_layers = opts$P1$nlayers,
      n_units = as.numeric((tune$results %>% select(starts_with("units_X_")))[1,]),
      penL1 = tune$results$penL1_X[1],
      penL2 = tune$results$penL2_X[1],
      lr = tune$results$lr_X[1]
    )
    autoencoder_Y <- statMatch.ACCA.autoencoder(
      n_var = ncol(data$mtx.don.NCV.raw),
      n_lat = tune$results$lat_Y[1],
      n_hidden_layers = opts$P1$nlayers,
      n_units = as.numeric((tune$results %>% select(starts_with("units_Y_")))[1,]),
      penL1 = tune$results$penL1_Y[1],
      penL2 = tune$results$penL2_Y[1],
      lr = tune$results$lr_Y[1]
    )
    fit_X <- autoencoder_X$model %>%
      fit(
        x = data$mtx.don.CV.raw,
        y = data$mtx.don.CV.raw,
        batch_size = opts$P1$batch_size,
        epochs = opts$P1$epochs,
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
        batch_size = opts$P1$batch_size,
        epochs = opts$P1$epochs,
        # validation_data = list(
        #   input$raw_data[[fold]]$mtx.don.CV.raw,
        #   input$raw_data[[fold]]$mtx.don.CV.raw,
        #   matrix(input$weights[[fold]])
        # ),
        verbose = 0,
        callbacks = opts$callbacks_keras,
        sample_weight = matrix(don.weights)
      )
    # Compute the CCA
    cca_raw <- statMatch.ACCA.phase1.cca(autoencoder_X, autoencoder_Y,
                                              data$mtx.rec.CV.raw,
                                              data$mtx.don.CV.raw,
                                              data$mtx.don.NCV.raw,
                                              don.weights,
                                              d=tune$results$d[1],
                                              opts)

    # Compute the predictions
    pred_raw <- statMatch.ACCA.phase1.predict(cca_raw,
                                              tune$results$h[1],
                                              data$mtx.don.NCV.raw,
                                              comp.mtx,
                                              don.weights,
                                              d=tune$results$d[1],
                                              opts)
    colnames(pred_raw) <- colnames(data$mtx.don.NCV.raw)
    pred <- mtx2df(pred_raw, var_names = names.NCV.cat)

    # Compute the sum of distances between the predictions and their compatible vectors in the donor data set
    #   This vector is used when averaging multiple models, for stability reasons.
    sum_dist_comp_donor <- compute_sum_dist_comp_donor(pred_raw, matrix(unlist(data$mtx.don.NCV.raw), nrow = nrow(comp.mtx)), comp.mtx)

    out <- list(tune = tune,
                fit_X = fit_X,
                fit_Y = fit_Y,
                df.match.don = pred,
                df.match = cbind(df.rec, pred),
                sum_dist_comp_donor = sum_dist_comp_donor)
  } else {
    if (length(names.NCV.cat) == 0)
      message("\nPhase 1 - No categorical variables to predict. Passing to the next phase ...")
    if (length(names.CV.cat) == 0 & length(names.NCV.cat) > 0)
      stop("No common categorical variables. Impossible to compute compatibility matrix.")
    out <- list(df.match = df.rec)
  }
  return(out)
}


compute_sum_dist_comp_donor <- function(pred, data_don, comp_mtx) {
  n_rec <- nrow(pred)
  out <- rep(NA, n_rec)
  for (i in seq_len(n_rec)) {
    out[i] <- sum(compute_eucl_dist_vec_mat_rcpp(pred[i,], data_don[comp_mtx[,i] == 1L, , drop = F]))
  }
  return(out)
}

