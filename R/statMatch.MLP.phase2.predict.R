#' @title statMatch.MLP.phase2.predict
#' @description Prediction in phase2 of the statistical matching algorithm using superOM
#' @param model model keras model, MLP trained using keras.
#' @param mtx.rec.CV.raw numeric matrix, common variables in the receiver data set (categorical variables transformed into dummies).
#' @param mtx.don.NCV.raw numeric matrix, non-common variables in the donor data set.
#' @param don.weights numeric vector, sample weights of the donor data set.
#' @param opts a list returned by the function \code{statMatch.MLP.options()} which contains all the options.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param keep_unconstrained a logical, whether or not the predictions are adjusted with the zero constraints.
#' @return A matrix, ...
#' @noRd

statMatch.MLP.phase2.predict <- function(model, mtx.rec.CV.raw, mtx.don.NCV.raw, don.weights, opts, zero.constraints, keep_unconstrained = F) {

  # Compute predictions from the fitted keras model
  PRED <- model %>%
    predict(
      x = mtx.rec.CV.raw,
      batch_size = opts$P2$batch_size,
      callbacks = opts$callbacks
    )

  # Go back to the original scale of the variables
  is.constrained.CV.rec <- startsWith(colnames(mtx.rec.CV.raw), "ZC_")
  if (opts$scaling == "z-score") {
    PRED <- scale_mtx_inv(PRED,
                          center = attr(mtx.don.NCV.raw, "scaled:center"),
                          scale = attr(mtx.don.NCV.raw, "scaled:scale"))
    if (length(zero.constraints) > 0) {
      ZC <- scale_mtx_inv(mtx.rec.CV.raw[, is.constrained.CV.rec],
                          center = attr(mtx.rec.CV.raw, "scaled:center")[is.constrained.CV.rec],
                          scale = attr(mtx.rec.CV.raw, "scaled:scale")[is.constrained.CV.rec])
    }
  } else if (opts$scaling == "min-max") {
    PRED <- normalize_mtx_inv(PRED,
                              min = attr(mtx.don.NCV.raw, "min"),
                              max = attr(mtx.don.NCV.raw, "max"))
    if (length(zero.constraints) > 0) {
      ZC <- normalize_mtx_inv(mtx.rec.CV.raw[, is.constrained.CV.rec],
                              min = attr(mtx.rec.CV.raw, "min")[is.constrained.CV.rec],
                              max = attr(mtx.rec.CV.raw, "max")[is.constrained.CV.rec])
    }
  }

  if (length(zero.constraints) > 0) ZC <- round(ZC) # <----------- Ensure that the zero constraints are equal 0 or 1 !!!

  names.NCV.don <- colnames(mtx.don.NCV.raw)
  colnames(PRED) <- names.NCV.don
  is.constrained <- names.NCV.don %in% zero.constraints
  if (!keep_unconstrained) {
    # Impose the constraint
    for (j in seq_len(length(names.NCV.don))) {
      if (is.constrained[j]) {
        constraint <- ZC[, paste0("ZC_", names.NCV.don[j], "...1")]
        PRED[constraint == 1, j] <- 0
      }
    }
  }

  return(PRED)
}
