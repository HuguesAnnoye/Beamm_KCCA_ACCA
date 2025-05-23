% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statMatch.ACCA.options.R
\name{statMatch.ACCA.options}
\alias{statMatch.ACCA.options}
\title{Statistical Matching: ACCA options}
\usage{
statMatch.ACCA.options(
  print.details = TRUE,
  scaling = c("z-score", "min-max", "no"),
  tuning_type = c("random"),
  d = NULL,
  p1_dmax = 2L,
  p1_dmin = 2L,
  p1_dstep = 2L,
  p2_dmax = 2L,
  p2_dmin = 2L,
  p2_dstep = 2L,
  p1_lat_X_min = 5L,
  p1_lat_X_max = 10L,
  p1_lat_X_step = 1L,
  p2_lat_X_min = 5L,
  p2_lat_X_max = 10L,
  p2_lat_X_step = 1L,
  p1_lat_Y_min = 5L,
  p1_lat_Y_max = 10L,
  p1_lat_Y_step = 1L,
  p2_lat_Y_min = 5L,
  p2_lat_Y_max = 10L,
  p2_lat_Y_step = 1L,
  p1_h = seq(0.1, 0.6, by = 0.05),
  p2_h = seq(0.1, 0.6, by = 0.05),
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
  p1_lr_min = -4,
  p1_lr_max = -1,
  p2_lr_min = -4,
  p2_lr_max = -1,
  p1_penL1_min = -6,
  p1_penL1_max = 0,
  p2_penL1_min = -6,
  p2_penL1_max = 0,
  p1_penL2_min = -6,
  p1_penL2_max = 0,
  p2_penL2_min = -6,
  p2_penL2_max = 0,
  p1_n_combs = 10L,
  p2_n_combs = 10L,
  n_fold = 5L,
  seed_phase1 = NULL,
  seed_phase2 = NULL,
  p1_man_params = NULL,
  p2_man_params = NULL,
  callbacks_keras = list(callback_model_checkpoint(filepath = paste0(tempfile(),
    ".hdf5"), monitor = "loss", verbose = 0, save_best_only = TRUE, save_weights_only =
    FALSE, mode = c("auto", "min", "max"), period = NULL, save_freq = "epoch"),
    callback_early_stopping(monitor = "loss", min_delta = 0, patience = 10, verbose = 0,
    restore_best_weights = TRUE))
)
}
\arguments{
\item{print.details}{a logical, whether or not details about the different procedures involved in the statistical matching are printed on the terminal.}

\item{scaling}{a character string, the name of the method used to scale the variables in the receiver and donor data sets: \code{"z-score"}, \code{"min-max"} or \code{"no"} (no scaling).}

\item{tuning_type}{a character string, the method used to find the optimal tuning parameters: \code{"random"} search only.}

\item{d}{a positive integer, number of latent variable used in CCA.}

\item{p1_dmax}{a positive integer, number of latent variable used in CCA.}

\item{p1_dmin}{a positive integer, number of latent variable used in CCA.}

\item{p1_dstep}{a positive integer, number of latent variable used in CCA.}

\item{p2_dmax}{a positive integer, number of latent variable used in CCA.}

\item{p2_dmin}{a positive integer, number of latent variable used in CCA.}

\item{p2_dstep}{a positive integer, number of latent variable used in CCA.}

\item{p1_lat_X_min}{a positive integer, minimum size of the latent space of the common variables in phase 1.}

\item{p1_lat_X_max}{a positive integer, maximum size of the latent space of the common variables in phase 1.}

\item{p1_lat_X_step}{a positive integer, step size of the latent space of the common variables in phase 1.}

\item{p2_lat_X_min}{a positive integer, minimum size of the latent space of the common variables in phase 2.}

\item{p2_lat_X_max}{a positive integer, maximum size of the latent space of the common variables in phase 2.}

\item{p2_lat_X_step}{a positive integer, step size of the latent space of the common variables in phase 2.}

\item{p1_lat_Y_min}{a positive integer, minimum size of the latent space of the non-common variables in phase 1.}

\item{p1_lat_Y_max}{a positive integer, maximum size of the latent space of the non-common variables in phase 1.}

\item{p1_lat_Y_step}{a positive integer, step size of the latent space of the non-common variables in phase 1.}

\item{p2_lat_Y_min}{a positive integer, minimum size of the latent space of the non-common variables in phase 2.}

\item{p2_lat_Y_max}{a positive integer, maximum size of the latent space of the non-common variables in phase 2.}

\item{p2_lat_Y_step}{a positive integer, step size of the latent space of the non-common variables in phase 2.}

\item{p1_h}{a numeric vector, possible values of the bandwidth parameter for the CCA procedure in phase 1.}

\item{p2_h}{a numeric vector, possible values of the bandwidth parameter for the CCA procedure in phase 2.}

\item{p1_nlayers}{a positive integer, number of hidden layers of the encoder and the decoder in phase 1.}

\item{p2_nlayers}{a positive integer, number of hidden layers of the encoder and the decoder in phase 2.}

\item{p1_epochs}{a positive integer, maximum number of epochs when training the autoencoder in phase 1.}

\item{p2_epochs}{a positive integer, maximum number of epochs when training the autoencoder in phase 2.}

\item{p1_batch_size}{a positive integer, number of samples per gradient update when training the autoencoder in phase 1.}

\item{p2_batch_size}{a positive integer, number of samples per gradient update when training the autoencoder in phase 2.}

\item{p1_u_min}{a positive integer, the minimum number of units in the hidden layers in phase 1.}

\item{p1_u_max}{a positive integer, the maximum number of units in the hidden layers in phase 1.}

\item{p1_u_step}{a positive integer, the step size for the number of units in the hidden layers in phase 1.}

\item{p2_u_min}{a positive integer, the minimum number of units in the hidden layers in phase 2.}

\item{p2_u_max}{a positive integer, the maximum number of units in the hidden layers in phase 2.}

\item{p2_u_step}{a positive integer, the step size for the number of units in the hidden layers in phase 2.}

\item{p1_lr_min}{a double, such that \code{10^-p1_lr_min} is the minimum value of the learning rate for the autoencoder in phase 1.}

\item{p1_lr_max}{a double, such that \code{10^-p1_lr_max} is the maximum value of the learning rate for the autoencoder in phase 1.}

\item{p2_lr_min}{a double, such that \code{10^-p2_lr_min} is the minimum value of the learning rate for the autoencoder in phase 2.}

\item{p2_lr_max}{a double, such that \code{10^-p2_lr_max} is the maximum value of the learning rate for the autoencoder in phase 2.}

\item{p1_penL1_min}{a double, such that \code{10^-p1_penL1_min} is the minimum value of the L1 penalty for the autoencoder in phase 1.}

\item{p1_penL1_max}{a double, such that \code{10^-p1_penL1_max} is the maximum value of the L1 penalty for the autoencoder in phase 1.}

\item{p2_penL1_min}{a double, such that \code{10^-p2_penL1_min} is the minimum value of the L1 penalty for the autoencoder in phase 2.}

\item{p2_penL1_max}{a double, such that \code{10^-p2_penL1_max} is the maximum value of the L1 penalty for the autoencoder in phase 2.}

\item{p1_penL2_min}{a double, such that \code{10^-p1_penL2_min} is the minimum value of the L2 penalty for the autoencoder in phase 1.}

\item{p1_penL2_max}{a double, such that \code{10^-p1_penL2_max} is the maximum value of the L2 penalty for the autoencoder in phase 1.}

\item{p2_penL2_min}{a double, such that \code{10^-p2_penL2_min} is the minimum value of the L2 penalty for the autoencoder in phase 2.}

\item{p2_penL2_max}{a double, such that \code{10^-p2_penL2_max} is the maximum value of the L2 penalty for the autoencoder in phase 2.}

\item{p1_n_combs}{a positive integer, number of random combinations of tuning parameters generated in phase 1 (used only when \code{tuning_type = "random"}).}

\item{p2_n_combs}{a positive integer, number of random combinations of tuning parameters generated in phase 2 (used only when \code{tuning_type = "random"}).}

\item{n_fold}{a positive integer, number of folds used for the cross-validation of the tuning parameters.}

\item{seed_phase1}{a positive integer, the seed for random number generation in phase 1 (a random value if \code{NULL}).}

\item{seed_phase2}{a positive integer, the seed for random number generation in phase 2 (a random value if \code{NULL}).}

\item{p1_man_params}{a data.frame, a manual way to specify the grid of tuning parameters in phase 1.}

\item{p2_man_params}{a data.frame, a manual way to specify the grid of tuning parameters in phase 2.}

\item{callbacks_keras}{a list containing the callbacks used by the keras model.}
}
\value{
A list with all the options required by the statistical matching algorithm based on ACCA.
}
\description{
Function to create the options for ACCA
}
\details{
DETAILS
}
