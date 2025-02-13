#' @title statMatch.ACCA.phase2.cca
#' @description CCA of phase 2 of the ACCA matching process
#' @param autoencoder_X keras model, autoencoder trained on the common variables.
#' @param autoencoder_Y keras model, autoencoder trained on the non-common variables.
#' @param mtx.rec.CV.raw numeric matrix, common variables in the receiver data set (categorical variables transformed into dummies).
#' @param mtx.don.CV.raw numeric matrix, common variables in the donor data set (categorical variables transformed into dummies).
#' @param mtx.don.NCV.raw numeric matrix, non-common variables in the donor data set (categorical variables transformed into dummies).
#' @param don.weights numeric vector, sample weights of the donor data set.
#' @param d numeric value, number of canonical variable to compute
#' @param opts a list returned by the function \code{statMatch.ACCA.options()} which contains all the options.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param keep_unconstrained a logical, whether or not the predictions are adjusted with the zero constraints.
#' @return RETURN
#' @details DETAILS
#' @importFrom magrittr '%>%'
#' @importFrom keras keras_model
#' @importFrom stats cancor
#' @importFrom pracma Rank
#' @noRd
#'
statMatch.ACCA.phase2.cca <- function(autoencoder_X, autoencoder_Y, mtx.rec.CV.raw, mtx.don.CV.raw,
                                          mtx.don.NCV.raw, don.weights, d, opts, zero.constraints, keep_unconstrained = F) {

  # Encode common variables (donor and receiver)
  model_X <- keras_model(inputs = autoencoder_X$inputs, outputs = autoencoder_X$encoder)
  pred_X_don_enc <- model_X %>% predict(x = mtx.don.CV.raw)
  pred_X_rec_enc <- model_X %>% predict(x = mtx.rec.CV.raw)

  # Encode non-common variables (donor)
  model_Y <- keras_model(inputs = autoencoder_Y$inputs, outputs = autoencoder_Y$encoder)
  pred_Y_don_enc <- model_Y %>% predict(x = mtx.don.NCV.raw)

  # Standardize the latent variables
  pred_Y_don_enc <- scale_mtx(pred_Y_don_enc, wt = don.weights)
  pred_X_don_enc <- scale_mtx(pred_X_don_enc, wt = don.weights)
  pred_X_rec_enc <- scale_mtx(
    pred_X_rec_enc,
    center = attr(pred_X_don_enc,"scaled:center"),
    scale = attr(pred_X_don_enc,"scaled:scale")
  )

  # Add column names (needed in case of perfect collinearity)
  colnames(pred_X_don_enc) <- paste0("V", 1:ncol(pred_X_don_enc))
  colnames(pred_X_rec_enc) <- paste0("V", 1:ncol(pred_X_rec_enc))

  # Normalize the donor weights
  weights_CCA <- sqrt(don.weights/sum(don.weights))
  # Fit CCA
  an.error.occured2 <- FALSE
  an.error.occured2 <<- FALSE

  fit_CCA <- tryCatch({cancor(weights_CCA * pred_X_don_enc, weights_CCA * pred_Y_don_enc, xcenter = FALSE, ycenter = FALSE)},
                      error = function(e) {an.error.occured2 <<- TRUE})

  # If the cancor fails then browser
  if (an.error.occured2) {
    if (Rank(pred_Y_don_enc) == 1) {
      can_var_CCA_mess <- "Perfect collinearity in all of the variables in the latent space of Y and non-calculable CCA."
      message(can_var_CCA_mess)
      can_var_CCA_don <- pred_X_don_enc
      can_var_CCA_rec <- pred_X_rec_enc
    } else if (Rank(pred_X_don_enc) == 1) {
      can_var_CCA_mess <- "Perfect collinearity in all of the variables in the latent space of X and non-calculable CCA."
      message(can_var_CCA_mess)
      can_var_CCA_don <- pred_X_don_enc
      can_var_CCA_rec <- pred_X_rec_enc
    } else if (Rank(pred_Y_don_enc) == 0) {
      can_var_CCA_mess <- "Null matrix in the latent space of Y and non-calculable CCA."
      message(can_var_CCA_mess)
      can_var_CCA_don <- pred_X_don_enc
      can_var_CCA_rec <- pred_X_rec_enc
    } else if (Rank(pred_X_don_enc) == 0) {
      can_var_CCA_mess <- "Null matrix in the latent space of X and non-calculable CCA."
      message(can_var_CCA_mess)
      can_var_CCA_don <- pred_X_don_enc
      can_var_CCA_rec <- pred_X_rec_enc
    } else {
      can_var_CCA_mess <- "Other unknown error and non-calculable CCA."
      message(can_var_CCA_mess)
      #browser()
      can_var_CCA_don <- pred_X_don_enc
      can_var_CCA_rec <- pred_X_rec_enc
    }
  } else {
    if (Rank(pred_Y_don_enc) == 1) {
      can_var_CCA_mess <- "Perfect collinearity in all of the variables in the latent space of Y but calculable CCA with reduced matrix Y."
      message(can_var_CCA_mess)
    }
  # Compute canonical variables
  if (nrow(fit_CCA$xcoef) == 1) {
    can_var_CCA_mess <- "Perfect collinearity in all of the variables in the latent space of X."
    message(can_var_CCA_mess)
    can_var_CCA_don <- pred_X_don_enc
    can_var_CCA_rec <- pred_X_rec_enc
  } else if (nrow(fit_CCA$xcoef) == ncol(pred_X_don_enc) | nrow(fit_CCA$xcoef) == ncol(pred_X_rec_enc)) {
    if (Rank(pred_Y_don_enc) == 1) {
      can_var_CCA_mess <- "Perfect collinearity in all of the variables in the latent space of Y but calculable CCA with reduced matrix Y."
    } else{
      can_var_CCA_mess <- "No-collinearity and calculabe CCA"
    }
    message(can_var_CCA_mess)
    pos_col <- min(d,ncol(fit_CCA$xcoef))
    can_var_CCA_don <- pred_X_don_enc %*% fit_CCA$xcoef[, 1:pos_col]
    can_var_CCA_rec <- pred_X_rec_enc %*% fit_CCA$xcoef[, 1:pos_col]
  } else {
    can_var_CCA_mess <- "Perfect collinearity in some of the variables in the latent space of X."
    message(can_var_CCA_mess)
    pos <- rownames(fit_CCA$xcoef)
    pos_col <- min(d,ncol(fit_CCA$xcoef))
    can_var_CCA_don <- pred_X_don_enc[,pos] %*% fit_CCA$xcoef[, 1:pos_col]
    can_var_CCA_rec <- pred_X_rec_enc[,pos] %*% fit_CCA$xcoef[, 1:pos_col]
  }
  }
  # Go back to the original scale of the variables
  if (opts$scaling == "z-score") {
    mtx.rec.CV.raw <- scale_mtx_inv(mtx.rec.CV.raw,
                                    center = attr(mtx.rec.CV.raw, "scaled:center"),
                                    scale = attr(mtx.rec.CV.raw, "scaled:scale"))
    mtx.don.CV.raw <- scale_mtx_inv(mtx.don.CV.raw,
                                    center = attr(mtx.don.CV.raw, "scaled:center"),
                                    scale = attr(mtx.don.CV.raw, "scaled:scale"))
    mtx.don.NCV.raw <- scale_mtx_inv(mtx.don.NCV.raw,
                                     center = attr(mtx.don.NCV.raw, "scaled:center"),
                                     scale = attr(mtx.don.NCV.raw, "scaled:scale"))
  } else if (opts$scaling == "min-max"  & length(zero.constraints) > 0) {
    mtx.rec.CV.raw <- normalize_mtx_inv(mtx.rec.CV.raw,
                                        min = attr(mtx.rec.CV.raw, "min"),
                                        max = attr(mtx.rec.CV.raw, "max"))
    mtx.don.CV.raw <- normalize_mtx_inv(mtx.don.CV.raw,
                                        min = attr(mtx.don.CV.raw, "min"),
                                        max = attr(mtx.don.CV.raw, "max"))
    mtx.don.NCV.raw <- normalize_mtx_inv(mtx.don.NCV.raw,
                                         min = attr(mtx.don.NCV.raw, "min"),
                                         max = attr(mtx.don.NCV.raw, "max"))
  }

  return(list(can_var_CCA_don=can_var_CCA_don,
              can_var_CCA_rec=can_var_CCA_rec,
              can_var_CCA_mess = can_var_CCA_mess,
              mtx.rec.CV.raw=mtx.rec.CV.raw,
              mtx.don.CV.raw=mtx.don.CV.raw,
              mtx.don.NCV.raw=mtx.don.NCV.raw))

}
