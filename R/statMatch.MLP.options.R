#' @title Statistical Matching: MLP options
#' @description Function to create the options for MLP
#'
#' @param print.details a logical, whether or not details about the different procedures involved in the statistical matching are printed on the terminal.
#' @param scaling a character string, the name of the method used to scale the variables in the receiver and donor data sets: \code{"z-score"}, \code{"min-max"} or \code{"no"} (no scaling).
#' @param tuning_type a character string, the method used to find the optimal tuning parameters: \code{"random"} search only.
#' @param p1_nlayers a positive integer, number of hidden layers of the MLP in phase 1.
#' @param p2_nlayers a positive integer, number of hidden layers of the MLP in phase 2.
#' @param p1_epochs a positive integer, maximum number of epochs when training the MLP in phase 1.
#' @param p2_epochs a positive integer, maximum number of epochs when training the MLP in phase 2.
#' @param p1_batch_size a positive integer, number of samples per gradient update when training the MLP in phase 1.
#' @param p2_batch_size a positive integer, number of samples per gradient update when training the MLP in phase 2.
#' @param p1_u_min a positive integer, the minimum number of units in the hidden layers in phase 1.
#' @param p1_u_max a positive integer, the maximum number of units in the hidden layers in phase 1.
#' @param p1_u_step a positive integer, the step size for the number of units in the hidden layers in phase 1.
#' @param p2_u_min a positive integer, the minimum number of units in the hidden layers in phase 2.
#' @param p2_u_max a positive integer, the maximum number of units in the hidden layers in phase 2.
#' @param p2_u_step a positive integer, the step size for the number of units in the hidden layers in phase 2.
#' @param p1_lr_min a double, such that \code{10^-p1_lr_min} is the minimum value of the learning rate for the MLP in phase 1.
#' @param p1_lr_max a double, such that \code{10^-p1_lr_max} is the maximum value of the learning rate for the MLP in phase 1.
#' @param p2_lr_min a double, such that \code{10^-p2_lr_min} is the minimum value of the learning rate for the MLP in phase 2.
#' @param p2_lr_max a double, such that \code{10^-p2_lr_max} is the maximum value of the learning rate for the MLP in phase 2.
#' @param p1_penL1_min a double, such that \code{10^-p1_penL1_min} is the minimum value of the L1 penalty for the MLP in phase 1.
#' @param p1_penL1_max a double, such that \code{10^-p1_penL1_max} is the maximum value of the L1 penalty for the MLP in phase 1.
#' @param p2_penL1_min a double, such that \code{10^-p2_penL1_min} is the minimum value of the L1 penalty for the MLP in phase 2.
#' @param p2_penL1_max a double, such that \code{10^-p2_penL1_max} is the maximum value of the L1 penalty for the MLP in phase 2.
#' @param p1_penL2_min a double, such that \code{10^-p1_penL2_min} is the minimum value of the L2 penalty for the MLP in phase 1.
#' @param p1_penL2_max a double, such that \code{10^-p1_penL2_max} is the maximum value of the L2 penalty for the MLP in phase 1.
#' @param p2_penL2_min a double, such that \code{10^-p2_penL2_min} is the minimum value of the L2 penalty for the MLP in phase 2.
#' @param p2_penL2_max a double, such that \code{10^-p2_penL2_max} is the maximum value of the L2 penalty for the MLP in phase 2.
#' @param p1_n_combs a positive integer, number of random combinations of tuning parameters generated in phase 1 (used only when \code{tuning_type = "random"}).
#' @param p2_n_combs a positive integer, number of random combinations of tuning parameters generated in phase 2 (used only when \code{tuning_type = "random"}).
#' @param n_fold a positive integer, number of folds used for the cross-validation of the tuning parameters.
#' @param seed_phase1 a positive integer, the seed for random number generation in phase 1 (a random value if \code{NULL}).
#' @param seed_phase2 a positive integer, the seed for random number generation in phase 2 (a random value if \code{NULL}).
#' @param p1_man_params a data.frame, a manual way to specify the grid of tuning parameters in phase 1.
#' @param p2_man_params a data.frame, a manual way to specify the grid of tuning parameters in phase 2.
#' @param callbacks_keras a list containing the callbacks used by the keras model.
#'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom stats runif
#' @importFrom keras callback_model_checkpoint callback_early_stopping
#'
#' @export

statMatch.MLP.options <- function(print.details = TRUE,
                                  scaling = c("z-score", "min-max", "no"),
                                  tuning_type = c("random"),
                                  p1_nlayers = 2L,
                                  p2_nlayers = 2L,
                                  p1_epochs = 200L,
                                  p2_epochs = 200L,
                                  p1_batch_size = 32L,
                                  p2_batch_size = 32L,
                                  p1_u_min = 20L,
                                  p1_u_max = 100L,
                                  p1_u_step = 5L,
                                  p2_u_min = 20L,
                                  p2_u_max = 100L,
                                  p2_u_step = 5L,
                                  p1_lr_min = -4,  # Attention in log-scale
                                  p1_lr_max = -1,  # Attention in log-scale
                                  p2_lr_min = -4,  # Attention in log-scale
                                  p2_lr_max = -1,  # Attention in log-scale
                                  p1_penL1_min = -6,  # Attention in log-scale
                                  p1_penL1_max = 0,  # Attention in log-scale
                                  p2_penL1_min = -6,  # Attention in log-scale
                                  p2_penL1_max = 0,  # Attention in log-scale
                                  p1_penL2_min = -6,  # Attention in log-scale
                                  p1_penL2_max = 0,  # Attention in log-scale
                                  p2_penL2_min = -6,  # Attention in log-scale
                                  p2_penL2_max = 0,  # Attention in log-scale
                                  # Argument for the tuning
                                  p1_n_combs = 10L,
                                  p2_n_combs = 10L,
                                  n_fold = 5L,
                                  seed_phase1 = NULL,
                                  seed_phase2 = NULL,
                                  p1_man_params = NULL,
                                  p2_man_params = NULL,
                                  callbacks_keras = list(
                                    # callback_reduce_lr_on_plateau(
                                    #   monitor = "val_loss",
                                    #   patience = 3
                                    # ),
                                    callback_model_checkpoint(
                                      filepath = paste0(tempfile(),".hdf5"),
                                      monitor = "loss",
                                      verbose = 0,
                                      save_best_only = TRUE,
                                      save_weights_only = FALSE,
                                      mode = c("auto", "min", "max"),
                                      period = NULL,
                                      save_freq = "epoch"
                                    ),
                                    callback_early_stopping(
                                      monitor = "loss",
                                      min_delta = 0,
                                      patience = 5,
                                      verbose = 0,
                                      restore_best_weights = TRUE
                                    )
                                  )) {
  scaling <- match.arg(scaling)
  tuning_type <- match.arg(tuning_type)
  if (length(print.details) != 1 | !is.logical(print.details) | any(is.na(print.details)))
    stop("MLP.options: argument 'print.details' should be a logical.")
  if (length(n_fold) != 1 | !is.integer(n_fold) | (n_fold <= 0))
    stop("ACCA: argument 'n_fold' should be a positive integer.")

  if (length(p1_nlayers) != 1 | !is.integer(p1_nlayers) | (p1_nlayers <= 0))
    stop("MLP.options: argument 'p1_nlayers' should be a positive integer.")
  if (length(p2_nlayers) != 1 | !is.integer(p2_nlayers) | (p2_nlayers <= 0))
    stop("MLP.options: argument 'p2_nlayers' should be a positive integer.")
  if (length(p1_epochs) != 1 | !is.integer(p1_epochs) | (p1_epochs <= 0))
    stop("MLP.options: argument 'p1_epochs' should be a positive integer.")
  if (length(p2_epochs) != 1 | !is.integer(p2_epochs) | (p2_epochs <= 0))
    stop("MLP.options: argument 'p2_epochs' should be a positive integer.")
  if (length(p1_batch_size) != 1 | !is.integer(p1_batch_size) | (p1_batch_size <= 0))
    stop("MLP.options: argument 'p1_batch_size' should be a positive integer.")
  if (length(p2_batch_size) != 1 | !is.integer(p2_batch_size) | (p2_batch_size <= 0))
    stop("MLP.options: argument 'p2_batch_size' should be a positive integer.")

  if (is.null(seed_phase1)) seed_phase1 <- sample(1:1e9, 1)
  else if (length(seed_phase1) != 1 | !is.integer(seed_phase1) | (seed_phase1 <= 0))
    stop("MLP.options: argument 'seed_phase1' should be a positive integer.")
  if (is.null(seed_phase2)) seed_phase2 <- sample(1:1e9, 1)
  else if (length(seed_phase2) != 1 | !is.integer(seed_phase2) | (seed_phase2 <= 0))
    stop("MLP.options: argument 'seed_phase2' should be a positive integer.")

  if (is.null(p1_man_params)) {
    if (length(p1_n_combs) != 1 | !is.integer(p1_n_combs) | (p1_n_combs <= 0))
      stop("MLP.options: argument 'p1_n_combs' should be a positive integer.")
    if (length(p1_u_min) != 1 | !is.integer(p1_u_min) | (p1_u_min <= 0))
      stop("MLP.options: argument 'p1_u_min' should be a positive integer.")
    if (length(p1_u_max) != 1 | !is.integer(p1_u_max) | (p1_u_max <= 0) | (p1_u_max < p1_u_min))
      stop("MLP.options: argument 'p1_u_max' should be a positive integer greater or equal to 'p1_u_min'.")
    if (length(p1_u_step) != 1 | !is.integer(p1_u_step) | (p1_u_step <= 0))
      stop("MLP.options: argument 'p1_u_step' should be a positive integer.")
    if (length(p1_lr_min) != 1 | !is.numeric(p1_lr_min) | (p1_lr_min > 0))
      stop("MLP.options: argument 'p1_lr_min' should be a double lower than zero.")
    if (length(p1_lr_max) != 1 | !is.numeric(p1_lr_max) | (p1_lr_max > 0) | (p1_lr_max < p1_lr_min))
      stop("MLP.options: argument 'p1_lr_max' should be a double lower than zero and greater or equal to 'p1_lr_min'.")
    if (length(p1_penL1_min) != 1 | !is.numeric(p1_penL1_min) | (p1_penL1_min > 0))
      stop("MLP.options: argument 'p1_penL1_min' should be a double lower than zero.")
    if (length(p1_penL1_max) != 1 | !is.numeric(p1_penL1_max) | (p1_penL1_max > 0) | (p1_penL1_max < p1_penL1_min))
      stop("MLP.options: argument 'p1_penL1_max' should be a double lower than zero and greater or equal to 'p1_penL1_min'.")
    if (length(p1_penL2_min) != 1 | !is.numeric(p1_penL2_min) | (p1_penL2_min > 0))
      stop("MLP.options: argument 'p1_penL2_min' should be a double lower than zero.")
    if (length(p1_penL2_max) != 1 | !is.numeric(p1_penL2_max) | (p1_penL2_max > 0) | (p1_penL2_max < p1_penL2_min))
      stop("MLP.options: argument 'p1_penL2_max' should be a double lower than zero and greater or equal to 'p1_penL2_min'.")
    if (tuning_type == "random") {
      hpar1 <- init_random_search_MLP(p1_n_combs, p1_lr_min, p1_lr_max, p1_penL1_min, p1_penL1_max,
                                      p1_penL2_min, p1_penL2_max, p1_nlayers, p1_u_min, p1_u_max, p1_u_step)
    }
  } else {
    hpar1 <- check_man_params_MLP(p1_man_params, p1_nlayers, "1")
  }

  if (is.null(p2_man_params)) {
    if (length(p2_n_combs) != 1 | !is.integer(p2_n_combs) | (p2_n_combs <= 0))
      stop("MLP.options: argument 'p2_n_combs' should be a positive integer.")
    if (length(p2_u_min) != 1 | !is.integer(p2_u_min) | (p2_u_min <= 0))
      stop("MLP.options: argument 'p2_u_min' should be a positive integer.")
    if (length(p2_u_max) != 1 | !is.integer(p2_u_max) | (p2_u_max <= 0) | (p2_u_max < p2_u_min))
      stop("MLP.options: argument 'p2_u_max' should be a positive integer greater or equal to 'p2_u_min'.")
    if (length(p2_u_step) != 1 | !is.integer(p2_u_step) | (p2_u_step <= 0))
      stop("MLP.options: argument 'p2_u_step' should be a positive integer.")
    if (length(p2_lr_min) != 1 | !is.numeric(p2_lr_min) | (p2_lr_min > 0))
      stop("MLP.options: argument 'p2_lr_min' should be a double lower than zero.")
    if (length(p2_lr_max) != 1 | !is.numeric(p2_lr_max) | (p2_lr_max > 0) | (p2_lr_max < p2_lr_min))
      stop("MLP.options: argument 'p2_lr_max' should be a double lower than zero and greater or equal to 'p2_lr_min'.")
    if (length(p2_penL1_min) != 1 | !is.numeric(p2_penL1_min) | (p2_penL1_min > 0))
      stop("MLP.options: argument 'p2_penL1_min' should be a double lower than zero.")
    if (length(p2_penL1_max) != 1 | !is.numeric(p2_penL1_max) | (p2_penL1_max > 0) | (p2_penL1_max < p2_penL1_min))
      stop("MLP.options: argument 'p2_penL1_max' should be a double lower than zero and greater or equal to 'p2_penL1_min'.")
    if (length(p2_penL2_min) != 1 | !is.numeric(p2_penL2_min) | (p2_penL2_min > 0))
      stop("MLP.options: argument 'p2_penL2_min' should be a double lower than zero.")
    if (length(p2_penL2_max) != 1 | !is.numeric(p2_penL2_max) | (p2_penL2_max > 0) | (p2_penL2_max < p2_penL2_min))
      stop("MLP.options: argument 'p2_penL2_max' should be a double lower than zero and greater or equal to 'p2_penL2_min'.")
    if (tuning_type == "random") {
      hpar2 <- init_random_search_MLP(p2_n_combs, p2_lr_min, p2_lr_max, p2_penL1_min, p2_penL1_max,
                                       p2_penL2_min, p2_penL2_max, p2_nlayers, p2_u_min, p2_u_max, p2_u_step)
    }
  } else {
    hpar2 <- check_man_params_MLP(p2_man_params, p2_nlayers, "2")
  }

  P1 <- list(
    nlayers = p1_nlayers,
    epochs = p1_epochs,
    batch_size = p1_batch_size,
    n_combs = p1_n_combs,
    hpar = hpar1
    # h_opt = p1_h
  )

  P2 <- list(
    nlayers = p2_nlayers,
    epochs = p2_epochs,
    batch_size = p2_batch_size,
    n_combs = p2_n_combs,
    hpar = hpar2
  )

  opts <- list(
    print.details = print.details,
    scaling = scaling,
    # Arguments for Keras
    P1 = P1,
    P2 = P2,
    # Argument for the tuning
    n_fold = n_fold,
    seed_phase1 = seed_phase1,
    seed_phase2 = seed_phase2,
    callbacks_keras = callbacks_keras
  )

  class(opts) <- "MLP.options"

  return(opts)
}



init_random_search_MLP <- function(n_combs, lr_min, lr_max, penL1_min, penL1_max, penL2_min, penL2_max,
                                   nlayers, u_min, u_max, u_step) {
  out <- data.frame(
    lr = 10^runif(n = n_combs, min = lr_min, max = lr_max),
    penL1 = 10^runif(n = n_combs, min = penL1_min, max = penL1_max),
    penL2 = 10^runif(n = n_combs, min = penL2_min, max = penL2_max)
  )
  for (i in seq_len(nlayers)) {
    out <- cbind(out, sample(seq(from = u_min, to = u_max, by = u_step), size = n_combs, replace = TRUE))
  }
  names(out) <- c("lr", "penL1", "penL2", paste0("units", seq_len(nlayers)))
  return(out)
}


check_man_params_MLP <- function(params, nlayers, phase = c("1", "2")) {
  phase <- match.arg(phase)
  if (!is.data.frame(params))
    stop("MLP.options: argument 'man_params_", phase, "' should be a data.frame object.")
  # Check if there is any missing column
  req_names <- c("lr", "penL1", "penL2", paste0("units", seq_len(nlayers)))
  check <- !req_names %in% colnames(params)
  missing <- paste0("'", req_names[check],"'", collapse = ", ")
  if (any(check))
    stop("MLP.options: column(s) ", missing, " missing in 'man_params_", phase, "' data.frame.")
  # Check numeric format
  check_names <- c("lr", "penL1", "penL2")
  for (i in 1:length(check_names)) {
    if (!is.numeric(params[,check_names[i]]) | any(params[,check_names[i]] <= 0))
      stop("MLP.options: column '", check_names[i], "' in 'man_params_", phase,
           "' data.frame should be a numeric vector with values grater than zero.")
  }
  # Check integer format
  check_names <- setdiff(req_names, check_names)
  for (i in 1:length(check_names)) {
    if (!is.integer(params[,check_names[i]]) | any(params[,check_names[i]] <= 0))
      stop("MLP.options: column '", check_names[i], "' in 'man_params_", phase,
           "' data.frame should be an integerO vector with values grater than zero.")
  }
  return(params)
}
