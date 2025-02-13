#' @title statMatch.KCCA.phase0
#' @description Statistical matching (phase0) using KCCA
#' @param df.rec a data.frame, the receiver data set with the varibales use in KCCA.
#' @param df.rec2 a data.frame, the receiver data set wit all the variables.
#' @param df.don a data.frame, the donor data set.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param opts a list returned by the function \code{statMatch.KCCA.one.phase.options()} which contains all the options.
#' @return RETURN
#' @details DETAILS
#' @importFrom dplyr select all_of starts_with
#' @importFrom magrittr '%>%'
#' @noRd


statMatch.KCCA.phase0 <- function(df.rec, df.don, df.rec2, don.weights, names.CV, names.NCV, zero.constraints, opts) {
  if (opts$print.details) {
    message("\nPhase 0 - Matching continuous variables.")
    message("\nCommon variables: ", paste0(names.CV, collapse = ", "), ".")
    message("\nNon-common variables: ", paste0(names.NCV, collapse = ", "), ".")
  }

  # Initialization
  names.CV.cat <- names.CV[unlist(lapply(df.don %>% select(all_of(names.CV)), is.factor))]
  names.NCV_cat <- names.NCV[unlist(lapply(df.don %>% select(all_of(names.NCV)), is.factor))]
  names.NCV_cont <- setdiff(names.NCV,names.NCV_cat) #names.NCV[unlist(lapply(df.don %>% select(all_of(names.NCV)), !is.factor))]
  if(length(c(names.NCV_cat,names.NCV_cont))!=length(names.NCV)) stop("names.NCV is not of the length of the sum of names.NCV_cat with names.NCV_cont.")
  # Tune the hyperparameters (if necessary)
  if ((opts$P0$n_combs) * (opts$P0$n_h) != 1) {
    if (isTRUE(opts$par)) {
      if (isTRUE(opts$par_fold)) {
      stop("Technique only exist in two phases")
      #tune <- statMatch.KCCA.phase0.tuning.par(df.don, names.CV, names.NCV, names.NCV_cat,names.NCV_cont, don.weights, zero.constraints, opts)
      } else {
      tune <- statMatch.KCCA.phase0.tuning.par2(df.don, names.CV, names.NCV, names.NCV_cat, names.NCV_cont, don.weights, zero.constraints, opts)
      }
    } else {
    tune <- statMatch.KCCA.phase0.tuning(df.don, names.CV, names.NCV, names.NCV_cat, names.NCV_cont, don.weights, zero.constraints, opts)
    }
  } else {
    tune <- list(results = opts$P0$hpar)
  }

  if (opts$print.details) {
    message("\nPhase 0. Fit and predict tuned model ...")
  }

  # Prepare data
  data <- list()
  data$mtx.don.CV.raw <- df2mtx(df.don %>% select(all_of(names.CV)))
  data$mtx.don.CV.raw_cont <- df2mtx(df.don %>% select(all_of(c(names.CV,names.NCV_cat))))
  data$mtx.don.NCV.raw <-  df2mtx(df.don %>% select(all_of(names.NCV)))
  data$mtx.don.NCV.raw_cat <- df2mtx(df.don %>% select(all_of(names.NCV_cat)))
  data$mtx.don.NCV.raw_cont <- df2mtx(df.don %>% select(all_of(names.NCV_cont)))
  data$mtx.rec.CV.raw <- df2mtx(df.rec %>% select(all_of(names.CV)))

  if (opts$scaling != "no") {
    data <- statMatch.scale.phase0(data, don.weights, opts$scaling)
  }

  # Fit and predict using the best hyperparameters

  fit <- statMatch.KCCA.fit(data$mtx.don.CV.raw,
    data$mtx.don.NCV.raw,
    weights = don.weights,
    X2 = data$mtx.rec.CV.raw,
    d = tune$results$d[1],
    h = tune$results$hx[1],
    hy = tune$results$hy[1],
    g = tune$results$g[1],
    rot = opts$rot
  )


  # Compute compatibility matrix
  comp.mtx <- compute_mtx_compatibilities(df.rec, df.don, names.CV.cat, opts$print.details, index = F)

  # Compute the predictions
  if (opts$type_predict == "matrix") {
    # Do the prediction using a matrix
    pred_raw <- statMatch.KCCA.phase1.predict(
      CV_X_A = fit$CV_X_A, # canonical variables for data set A
      CV_X_B = fit$CV_X_B, # canonical variables for data set B
      Y = data$mtx.don.NCV.raw_cat, # Non-common variables for data set A
      h = tune$results$h_cat[1],
      d = tune$results$d[1],
      kernel_predict = opts$P0$kernel_predict_cat,
      scaling = opts$scaling,
      rot = opts$rot,
      names.NCV = colnames(data$mtx.don.NCV.raw_cat),
      weights = don.weights,
      matrix.tot.possibilitities = comp.mtx,
      print.details = opts$print.details
    )
    pred_cat<- mtx2df(pred_raw, var_names = names.NCV_cat)

  } else if (opts$type_predict == "loop") {
    pred_raw <- statMatch.CCA.phase1.predict(
      CV_X_A = fit$CV_X_A, # canonical variable for data set A
      CV_X_B = fit$CV_X_B, # canonical variable for data set B
      h = tune$results$h_cat[1],
      Y = data$mtx.don.NCV.raw,
      comp.mtx = comp.mtx,
      don.weights = don.weights,
      opts = opts
    )
    colnames(pred_raw) <- colnames(data$mtx.don.NCV.raw)
    pred_cat<- mtx2df(pred_raw, var_names = names.NCV_cat)
  }

  #data$mtx.don.CV.raw <- df2mtx(cbind(df.don %>% select(all_of(c(names.CV,names.NCV_cat)))))
  data$mtx.rec.CV.raw_cont <- df2mtx(cbind(df.rec, pred_cat))
  if (opts$scaling == "z-score") {
  data$mtx.rec.CV.raw_cont <- scale_mtx(data$mtx.rec.CV.raw_cont,
                                    center = attr(data$mtx.don.CV.raw_cont, "scaled:center"),
                                    scale = attr(data$mtx.don.CV.raw_cont, "scaled:scale"))
  } else if (opts$scaling ==  "min-max"){
    data$mtx.rec.CV.raw_cont <- normalize_mtx(data$mtx.rec.CV.raw_cont,
                                          min = attr(data$mtx.don.CV.raw_cont, "min"),
                                          max = attr(data$mtx.don.CV.raw_cont, "max"))
  }

  # Compute the predictions
  if (opts$type_predict == "matrix") {
    pred_raw <- statMatch.KCCA.phase2.predict(
      CV_X_A = fit$CV_X_A, # canonical variables for data set A
      CV_X_B = fit$CV_X_B, # canonical variables for data set B
      Y = data$mtx.don.NCV.raw_cont, # Non-common variables for data set A
      h = tune$results$h[1],
      d = tune$results$d[1],
      kernel_predict = opts$P0$kernel_predict,
      scaling = opts$scaling,
      rot = opts$rot,
      names.NCV = colnames(data$mtx.don.NCV.raw_cont),
      weights = don.weights,
      data_cat_rec = cbind(df.rec, pred_cat),
      data_cat_don = df.don,
      zero.constraints = zero.constraints,
      matrix.tot.possibilitities = NULL,
      print.details = opts$print.details
    )
  } else if (opts$type_predict == "loop") {
    # Do the prediction using a rcpp loop
    pred_raw <- statMatch.CCA.phase2.predict(
      CV_X_A = fit$CV_X_A, # canonical variables for data set A
      CV_X_B = fit$CV_X_B, # canonical variables for data set B
      h = tune$results$h[1],
      mtx.rec.CV.raw = data$mtx.rec.CV.raw,
      mtx.don.CV.raw = data$mtx.don.CV.raw_cont,
      mtx.don.NCV.raw = data$mtx.don.NCV.raw_cont,
      don.weights = don.weights,
      opts = opts,
      zero.constraints = zero.constraints,
      keep_unconstrained = T
    )
    colnames(pred_raw) <- colnames(data$mtx.don.NCV.raw_cont)
  }
  pred_cont<- mtx2df(pred_raw, var_names = names.NCV_cont)
  out <- list(
    tune = tune,
    fit = fit,
    df.match.don = cbind(pred_cat, pred_cont),
    df.match = cbind(df.rec2, pred_cat, pred_cont)
  )
  return(out)
}
