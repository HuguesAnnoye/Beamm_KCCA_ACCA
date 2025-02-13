#' @title MULTILAYER PERCEPTRON PHASE 1 PREDICTION
#' @description Prediction of phase 1 of the MLP matching process
#' @param model keras model, MLP trained using keras.
#' @param mtx.rec.CV.raw numeric matrix, common variables in the receiver data set (categorical variables transformed into dummies).
#' @param mtx.don.NCV.raw numeric matrix, non-common variables in the donor data set (categorical variables transformed into dummies).
#' @param comp.mtx integer matrix of size N.don-by-N.rec containing a value 1L if an individual in the donor data set is compatible with another one in the receiver data set.
#' @param don.weights numeric vector, sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.MLP.options()} which contains all the options.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @importFrom magrittr '%>%'
#' @noRd

statMatch.MLP.phase1.predict <- function(model, mtx.rec.CV.raw, mtx.don.NCV.raw, comp.mtx, don.weights, opts) {

  # Compute predictions from the keras model
  pred_Y <- model %>%
    predict(
      x = mtx.rec.CV.raw,
      batch_size = opts$P1$batch_size,
      callbacks = opts$callbacks
    )
  N.rec <- ncol(comp.mtx)
  N.don <- nrow(comp.mtx)
  # Transform keras predictions and donor NCV from list to matrix
  pred_Y <- matrix(unlist(pred_Y), nrow = N.rec)
  mtx.don.NCV.raw <- matrix(unlist(mtx.don.NCV.raw), nrow = N.don)
  # Based on the predictions, pick up at random among the compatible individuals in the donor
  idx <- 1:N.don
  mtx_PRED <- matrix(NA_real_, N.rec, ncol(pred_Y))
  for (i in seq_len(N.rec)) {
    prob_i <- don.weights
    idx.don.comp <- idx[comp.mtx[,i] == 1L]
    # Compute the euclidean distance to the compatible individuals in the donor data set
    distances <- compute_eucl_dist_vec_mat_rcpp(pred_Y[i,], mtx.don.NCV.raw[idx.don.comp, , drop = F])
    # Compute the gaussian kernel
    N.tmp <- length(idx.don.comp)
    if (N.tmp > 1)
      sigma <- wtd_sd_rcpp(distances, don.weights[idx.don.comp])
    else sigma <- 1
    h <- 1.06 * sigma * N.tmp ^ (-1/5)
    h_min <- min(distances)
    h_final <- pmax(h, h_min)
    if (h_final == 0) h_final <- 1
    kern_dist <-  sqrt(2*pi) * exp(-((distances/h_final)^2)/2)
    # Final weights used to pick up at random
    prob_i[idx.don.comp] <- kern_dist * prob_i[idx.don.comp]
    prob_i[-idx.don.comp] <- 0
    # Pick up an individual at random
    mtx_PRED[i,] <- mtx.don.NCV.raw[sample(N.don, 1, prob = prob_i),]
  }
  return(mtx_PRED)
}
