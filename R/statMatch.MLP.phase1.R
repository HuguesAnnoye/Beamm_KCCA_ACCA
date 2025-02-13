#' @title statMatch.MLP.phase1
#' @description Statistical matching (phase1) using MMLP.
#' @param df.rec a data.frame, the recipient data set.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param opts a list, supplied using the function \code{statMatch.superOM.options()}, containing all the options.
#' @return A list, ...
#'
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom keras k_clear_session
#' @importFrom tensorflow set_random_seed
#' @noRd

statMatch.MLP.phase1 <- function(df.rec, df.don, don.weights, names.CV, names.NCV, opts) {

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
      tune <- statMatch.MLP.phase1.tuning(df.don, names.CV, names.NCV.cat, don.weights, opts)
    else tune <- list(results = opts$P1$hpar)

    if (opts$print.details)
      message("\nPhase 1. Fit and predict tuned model ...")

    # Prepare data
    data <- list()
    data$mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)))
    data$mtx.don.NCV.raw_tmp <- df2mtx(df.don %>% select(all_of(names.NCV.cat)))
    data$mtx.rec.CV.raw <- df2mtx(df.rec %>% select(all_of(names.CV)))
    names.NCV.raw <- colnames(data$mtx.don.NCV.raw_tmp)
    data$mtx.don.NCV.raw <- list()
    for (i in seq_len(length(names.NCV.cat))) {
      sel <- startsWith(names.NCV.raw, names.NCV.cat[i])
      data$mtx.don.NCV.raw[[i]] <- data$mtx.don.NCV.raw_tmp[, sel, drop = F]
    }
    data$mtx.don.NCV.raw_tmp <- NULL

    if (opts$scaling != "no")
      data <- statMatch.scale(data, don.weights, opts$scaling, scale_NCV = F)

    # Compute compatibility matrix
    comp.mtx <- compute_mtx_compatibilities(df.rec, df.don, names.CV.cat, opts$print.details, index = F)

    # Fit and predict using the best hyperparameters
    # Clear previous keras session
    k_clear_session()
    # Set same seed used in cross-validation
    set_random_seed(opts$seed_phase1)
    # Initialize the MMLP with keras
    model_phase1 <- statMatch.MLP.phase1.model(
      n_X = ncol(data$mtx.don.CV.raw),
      n_Y = length(data$mtx.don.NCV.raw),
      n_levels_y = sapply(data$mtx.don.NCV.raw, ncol),
      n_hidden_layers = opts$P1$nlayers,
      n_units = as.numeric((tune$results %>% select(starts_with("units")))[1,]),
      penL1 = tune$results$penL1[1],
      penL2 = tune$results$penL2[1],
      learning_rate = tune$results$lr[1]
    )
    # Fit model with keras
    fit <- model_phase1 %>%
      fit(
        x = data$mtx.don.CV.raw,
        y = data$mtx.don.NCV.raw,
        batch_size = opts$P1$batch_size,
        epochs = opts$P1$epochs,
        # validation_data = list(
        #   data$mtx.don.CV.raw,
        #   data$mtx.don.NCV.raw,
        #   matrix(don.weights)
        # ),
        verbose = 0,
        callbacks = opts$callbacks_keras,
        sample_weight = matrix(don.weights)
      )
    pred_raw <- statMatch.MLP.phase1.predict(model_phase1,
                                             data$mtx.rec.CV.raw,
                                             data$mtx.don.NCV.raw,
                                             comp.mtx,
                                             don.weights,
                                             opts)
    colnames(pred_raw) <- names.NCV.raw
    pred <- mtx2df(pred_raw, var_names = names.NCV.cat)

    # Compute the sum of distances between the predictions and their compatible vectors in the donor data set
    #   This vector is used when averaging multiple models, for stability reasons.
    sum_dist_comp_donor <- compute_sum_dist_comp_donor(pred_raw, matrix(unlist(data$mtx.don.NCV.raw), nrow = nrow(comp.mtx)), comp.mtx)

    out <- list(tune = tune,
                fit = fit,
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
