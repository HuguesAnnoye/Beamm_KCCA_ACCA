#' @title statMatch.superOM.phase2.predict
#' @description Prediction in phase2 of the statistical matching algorithm using superOM
#' @param fit an object of class "kohonen".
#' @param data a list containing three data.frames: \code{mtx.don.CV.raw}, \code{mtx.don.NCV.raw} and \code{mtx.rec.CV.raw}.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param scaling a character string, the name of the method used to scale the variables in the receiver and donor data sets: \code{"z-score"}, \code{"min-max"} or \code{"no"} (no scaling).
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param keep_unconstrained a logical, whether or not the predictions are adjusted with the zero constraints.
#' @return A matrix, ...
#'
#' @importFrom stats predict
#' @noRd

statMatch.superOM.phase2.predict <- function(fit, data, don.weights, scaling, zero.constraints, keep_unconstrained = F) {

  N.rec <- nrow(data$mtx.rec.CV.raw)
  # N.don <- nrow(data$mtx.don.CV.raw)
  N.NCV.don <- ncol(data$mtx.don.NCV.raw)
  # N.CV.don <- ncol(data$mtx.don.CV.raw)

  # Predict kohonen
  list.NEW.raw <- list(layer1 = data$mtx.rec.CV.raw,
                       layer2 = matrix(NA_real_, N.rec, N.NCV.don))
  pred.superOM <- predict(fit, newdata = list.NEW.raw, whatmap = "layer1")

  # Go back to the original scale of the variables
  if (scaling == "z-score") {
    mtx.don.CV.raw <- scale_mtx_inv(data$mtx.don.CV.raw,
                                    center = attr(data$mtx.don.CV.raw, "scaled:center"),
                                    scale = attr(data$mtx.don.CV.raw, "scaled:scale"))
    mtx.don.NCV.raw <- scale_mtx_inv(data$mtx.don.NCV.raw,
                                     center = attr(data$mtx.don.NCV.raw, "scaled:center"),
                                     scale = attr(data$mtx.don.NCV.raw, "scaled:scale"))
    mtx.rec.CV.raw <- scale_mtx_inv(data$mtx.rec.CV.raw,
                                    center = attr(data$mtx.rec.CV.raw, "scaled:center"),
                                    scale = attr(data$mtx.rec.CV.raw, "scaled:scale"))
  } else if (scaling == "min-max") {
    mtx.don.CV.raw <- normalize_mtx_inv(data$mtx.don.CV.raw,
                                        min = attr(data$mtx.don.CV.raw, "min"),
                                        max = attr(data$mtx.don.CV.raw, "max"))
    mtx.don.NCV.raw <- normalize_mtx_inv(data$mtx.don.NCV.raw,
                                         min = attr(data$mtx.don.NCV.raw, "min"),
                                         max = attr(data$mtx.don.NCV.raw, "max"))
    mtx.rec.CV.raw <- normalize_mtx_inv(data$mtx.rec.CV.raw,
                                        min = attr(data$mtx.rec.CV.raw, "min"),
                                        max = attr(data$mtx.rec.CV.raw, "max"))
  }

  is.constrained.CV.don <- startsWith(colnames(mtx.don.CV.raw), "ZC_")                      #
  is.constrained.CV.rec <- startsWith(colnames(mtx.rec.CV.raw), "ZC_")                      #
  mtx.don.CV.raw[, is.constrained.CV.don] <- round(mtx.don.CV.raw[, is.constrained.CV.don]) # <------------------ Ensure that the zero constraints are equal 0 or 1 !!!
  mtx.rec.CV.raw[, is.constrained.CV.rec] <- round(mtx.rec.CV.raw[, is.constrained.CV.rec]) #

  n.neurons <- fit$grid$xdim * fit$grid$ydim
  pred.cb.neuron <- matrix(NA, n.neurons, N.NCV.don)
  names.NCV.don <- colnames(data$mtx.don.NCV.raw)
  is.constrained <- names.NCV.don %in% zero.constraints
  for (neuron in 1:n.neurons) {
    idx <- fit$unit.classif == neuron
    for (j in 1:N.NCV.don) {
      wt <- don.weights
      if (is.constrained[j]) {
        constraint <- mtx.don.CV.raw[, paste0("ZC_", names.NCV.don[j], "...1")]
        wt[!idx | constraint == 1] <- 0
      } else {
        wt[!idx] <- 0
      }
      if (all(wt == 0)) {
        pred.cb.neuron[neuron, j] <- 0
      } else {
        pred.cb.neuron[neuron, j] <- wtd_mean_rcpp(mtx.don.NCV.raw[, j], wt)
      }
    }
  }

  PRED <- pred.cb.neuron[pred.superOM$unit.classif, ]
  colnames(PRED) <- names.NCV.don
  if (!keep_unconstrained) {
    # Impose the constraint
    for (j in 1:N.NCV.don) {
      if (is.constrained[j]) {
        constraint <- mtx.rec.CV.raw[, paste0("ZC_", names.NCV.don[j], "...1")]
        PRED[constraint == 1, j] <- 0
      }
    }
  }

  return(PRED)
}
