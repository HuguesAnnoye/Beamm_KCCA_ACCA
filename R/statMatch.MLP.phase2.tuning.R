#' @title statMatch.MLP.phase2.tuning
#' @description Tuning (phase2) of the statistical matching algorithm using MLP
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param opts a list returned by the function \code{statMatch.ACCA.options()} which contains all the options.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom magrittr %>%
#' @importFrom keras k_clear_session
#' @importFrom tensorflow set_random_seed
#' @noRd

statMatch.MLP.phase2.tuning <- function(df.don, names.CV, names.NCV, don.weights, zero.constraints, opts) {

  # Prepare data for cross-validation
  input <- prepare_data_tuning(df.don, names.CV, names.NCV, don.weights, opts)
  n_X <- ncol(input$raw_data[[1]]$mtx.don.CV.raw)
  n_Y <- length(names.NCV)

  results <- opts$P2$hpar
  results$val_objfun <- NA
  results$avg_epochs <- NA
  results$time <- NA
  info <- list(epochs = matrix(NA, opts$P2$n_combs, opts$n_fold))

  message("\nPhase 2. Tuning MLP ...")

  # Perform cross-validation
  for (i in seq_len(opts$P2$n_combs)) {
    t_start <- Sys.time()
    # Clear previous keras session
    k_clear_session()
    # Set the same seed for each combination of hyperparameters
    # set.seed(opts$seed_phase2)
    set_random_seed(opts$seed_phase2)
    # Initialize the MMLP with keras
    model <- statMatch.MLP.phase2.model(
      n_X, n_Y,
      n_hidden_layers = opts$P2$nlayers,
      n_units = as.numeric((opts$P2$hpar %>% select(starts_with("units")))[i,]),
      penL1 = opts$P2$hpar$penL1[i],
      penL2 = opts$P2$hpar$penL2[i],
      learning_rate = opts$P2$hpar$lr[i]
    )

    # Predict and stack results over the k folds
    PRED.raw <- NULL
    for (fold in seq_len(opts$n_fold)) {
      message(Sys.time(), ". Combination ", i, " of ", opts$P2$n_combs, ". Fold ", fold, " of ", opts$n_fold)
      # Fit model with keras
      history <- model %>%
        fit(
          x = input$raw_data[[fold]]$mtx.don.CV.raw,
          y = input$raw_data[[fold]]$mtx.don.NCV.raw,
          batch_size = opts$P2$batch_size,
          epochs = opts$P2$epochs,
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
      pred <- statMatch.MLP.phase2.predict(model,
                                           input$raw_data[[fold]]$mtx.rec.CV.raw,
                                           input$raw_data[[fold]]$mtx.don.NCV.raw,
                                           input$weights[[fold]],
                                           opts, zero.constraints)
      # Stack the predictions
      PRED.raw <- rbind(PRED.raw, pred)

    }
    # Store elapsed time
    results$time[i] <-  as.numeric(difftime(Sys.time(), t_start, units = 'mins'))
    # Compute average of the epochs over the k folds
    results$avg_epochs[i] <- mean(info$epochs[i,])
    # Compute the objective function
    results$val_objfun[i] <- compute_obj_function(input$mtx.don.NCV.raw, PRED.raw, don.weights, type = "wsRMSE")
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


# statMatch.MLP.phase2.tuning <- function(df.don, df.rec, data, weights = NULL, job = NULL,
#                                         zero.constraints = NULL, callback = NULL,
#                                         names = NULL, opts = NULL) {
#
#   if (opts$print.details) {message("Part 2: Fitting")}
#
#   kfold <- createFolds(df.don[,1], k = opts$n_fold)
#
#   cv_mse_cont <- matrix(NA, nrow = nrow(opts$P2$hpar), ncol = 1)
#
#   for (j in 1:nrow(opts$P2$hpar)){
#
#     pred <- NULL
#
#     MLP <- keras_model_sequential() %>%
#           layer_dense(input_shape = ncol(data$don.CV), units = opts$P2$hpar[j, 1], activation = 'relu', activity_regularizer = regularizer_l1(l = opts$P2$hpar[j, 4])) %>%
#           layer_batch_normalization() %>%
#           layer_dropout(rate = 0.2) %>%
#           layer_dense(units = opts$P2$hpar[j, 2], activation = 'relu', activity_regularizer = regularizer_l1(l = opts$P2$hpar[j, 4])) %>%
#           layer_dropout(rate = 0.1) %>%
#           layer_dense(units = ncol(data$don.NCV), activation = 'relu') %>%
#           compile(loss = opts$P2$loss, weighted_metrics = opts$P2$metric,
#                   optimizer = optimizer_adam(lr = opts$P2$hpar[j, 3]))
#
#     for (i in 1:length(kfold)){
#
#       cv_Wps <- as.matrix(weights$Wps) %>% .[-kfold[[i]],]
#       cv_Wpst <- as.matrix(weights$Wps) %>% .[kfold[[i]],]
#
#       X <- data$don.CV[-kfold[[i]],]
#       Xt <- data$don.CV[kfold[[i]],]
#       Y <- data$don.NCV[-kfold[[i]],]
#       Xot <- data$don.CV[kfold[[i]],]
#       Yot <- data$don.NCV[kfold[[i]],]
#
#       history <- MLP %>%
#             keras::fit(X, Y, epochs = 10, verbose = 0, validation_split = 0.2,
#                        callbacks = callback, sample_weight = cv_Wps)
#
#       Yt_hat <- MLP %>%
#             keras::predict_on_batch(x = Xt)
#
#       colnames(Yt_hat) <- colnames(Y)
#
#       tmp <- Yot
#
#       pred.tot <- as.data.frame(cbind(Xot, Yt_hat))
#       for (name in zero.constraints){
#         ifelse(pred.tot[, paste0("ZC_", name, "...1")] == 1, 0, name)
#       }
#       mtx_PRED <- as.matrix(pred.tot[, names$NCV_cont])
#       pred <- rbind(pred, mtx_PRED)
#       print(paste("i = ", i))
#     }
#     cv_mse_cont[j] <- compute_obj_function(data$don.NCV, pred, weights$Wps, type = opts$objmethod_cont)
#     print(paste("j = ", j))
#   }
#
#   mini <- which(cv_mse_cont == min(cv_mse_cont), arr.ind = TRUE)
#
#   return(mini)
# }
