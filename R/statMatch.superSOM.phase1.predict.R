#' @title statMatch.superOM.phase1.predict
#' @description Prediction in phase1 of the statistical matching algorithm using superOM
#' @param fit an object of class "kohonen".
#' @param data a list containing three data.frames: \code{mtx.don.CV.raw}, \code{mtx.don.NCV.raw} and \code{mtx.rec.CV.raw}.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param comp.idx, a matrix with 2 columns, the positions of the compatible individuals between the donor and receiver data sets.
#' @param scaling a character string, the name of the method used to scale the variables in the receiver and donor data sets: \code{"z-score"}, \code{"min-max"} or \code{"no"} (no scaling).
#' @return A matrix, ...
#'
#' @importFrom stats predict
#' @noRd

statMatch.superOM.phase1.predict <- function(fit, data, don.weights, comp.idx, scaling) {
  N.rec <- nrow(data$mtx.rec.CV.raw)
  N.don <- nrow(data$mtx.don.CV.raw)
  N.NCV.don <- ncol(data$mtx.don.NCV.raw)

  list.NEW.raw <- list(layer1 = data$mtx.rec.CV.raw,
                       layer2 = matrix(NA_real_, N.rec, N.NCV.don))
  pred.superOM <- predict(fit, newdata = list.NEW.raw, whatmap = "layer1")

  idx <- 1:N.don
  PRED_raw <- matrix(NA_real_, N.rec, N.NCV.don)
  for (i in 1:N.rec) {
    prob_i <- don.weights
    idx.don.neuron <- which(pred.superOM$unit.classif[i] == fit$unit.classif)
    idx.don.comp <- idx[comp.idx[,i] == 1L]
    idx.don.neuron.comp <- intersect(idx.don.neuron, idx.don.comp)
    # Check if there are compatible individuals in the predicted neuron
    if (length(idx.don.neuron.comp) > 0) {
      # Compute the euclidean distance to the compatible individuals in the donor data set:
      distances <- compute_eucl_dist_vec_mat_rcpp(data$mtx.rec.CV.raw[i,], data$mtx.don.CV.raw[idx.don.neuron.comp, , drop = F])
      # Compute the gaussian kernel
      N.tmp <- length(idx.don.neuron.comp)
      if (N.tmp > 1)
        sigma <- wtd_sd_rcpp(distances, don.weights[idx.don.neuron.comp])
      else sigma <- 1
      h <- 1.06 * sigma * N.tmp ^ (-1/5)
      h_min <- min(distances)
      h_final <- pmax(h, h_min)
      if (h_final == 0) h_final <- 1
      kern_dist <-  sqrt(2*pi) * exp(-((distances/h_final)^2)/2)
      prob_i[idx.don.neuron.comp] <- kern_dist * prob_i[idx.don.neuron.comp]
      prob_i[-idx.don.neuron.comp] <- 0
    } else {
      # In this case there are no compatible individuals in the predicted neuron.
      # Compute the euclidean distance to the compatible individuals in the donor data set:
      distances <- compute_eucl_dist_vec_mat_rcpp(data$mtx.rec.CV.raw[i,], data$mtx.don.CV.raw[idx.don.comp, , drop = F])
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
      prob_i[idx.don.comp] <- kern_dist * prob_i[idx.don.comp]
      prob_i[-idx.don.comp] <- 0
    }
    PRED_raw[i,] <- data$mtx.don.NCV.raw[sample(N.don, 1, prob = prob_i),]
  }

  if (scaling == "z-score") {
    PRED_raw <- scale_mtx_inv(PRED_raw,
                              center = attr(data$mtx.don.NCV.raw,"scaled:center"),
                              scale = attr(data$mtx.don.NCV.raw,"scaled:scale"))
  } else if (scaling == "min-max") {
    PRED_raw <- normalize_mtx_inv(PRED_raw,
                                  min = attr(data$mtx.don.NCV.raw,"min"),
                                  max = attr(data$mtx.don.NCV.raw,"max"))
  }

  PRED_raw <- round(PRED_raw) # <------------------------------------- Ensure that the matrix cantains only 0s or 1s !!!

  colnames(PRED_raw) <- colnames(data$mtx.don.NCV.raw)

  return(PRED_raw)
}
