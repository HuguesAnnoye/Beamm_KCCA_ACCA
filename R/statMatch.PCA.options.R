#' @title Statistical Matching: PCA options
#' @description Function to create the options for PCA
#'
#' @param print.details a logical, print details or not (in the tuning part only print the value of the objective function for each parameters combination)
#' @param print.details_tuning a logical, if true print all the details in the tuning part
#' @param tuning_type a character string, the method used to find the optimal tuning parameters: \code{"random"}, \code{"two h"} or \code{"one h"}
#' @param scaling a character string, the name of the method used to scale the variables in the receiver and donor data sets: \code{"z-score"}, \code{"min-max"} or \code{"no"} (no scaling).
#' @param d a positive integer, number of latent variable used in CCA.
#' @param p1_dmax a positive integer, number of latent variable used in CCA.
#' @param p1_dmin a positive integer, number of latent variable used in CCA.
#' @param p1_dstep a positive integer, number of latent variable used in CCA.
#' @param p2_dmax a positive integer, number of latent variable used in CCA.
#' @param p2_dmin a positive integer, number of latent variable used in CCA.
#' @param p2_dstep a positive integer, number of latent variable used in CCA.
#' @param rot a logical, if TRUE the bandwidth is h multiply by the variance.
#' @param type_predict a character string, Type of prediction \code{"matrix"} or \code{"loop"}.
#' @param n_fold a positive integer, number of folds used for the cross-validation of the tuning parameters.
#' @param names.CV_Cat_prio a character vector, Name of common categorical variable that have to match.
#' @param p1_kernel_predict a character string, Type of kernel use for prediction for categorical variable.
#' @param p1_kernel_predict_tuning a character string, Type of kernel use for the tuning for categorical variable by default the same than kernel_predict.
#' @param p2_kernel_predict a character string, Type of kernel use for prediction for continuous variable.
#' @param p1_h a numeric vector, the vector of value for h, the bandwidth of the Matching kernel for categorical variable.
#' @param p1_hmax a numeric value, the maximal value for h, the bandwidth of the Matching kernel for categorical variable.
#' @param p1_hmin a numeric value, the minimal value for h, the bandwidth of the Matching kernel for categorical variable.
#' @param p1_hstep a numeric value, a step for the grid of h for categorical variable.
#' @param p2_h a numeric vector, the vector of value for h, the bandwidth of the Matching kernel for continuous variable.
#' @param p2_hmax a numeric value, the maximal value for h, the bandwidth of the Matching kernel for continuous variable.
#' @param p2_hmin a numeric value, the minimal value for h, the bandwidth of the Matching kernel for continuous variable.
#' @param p2_hstep a numeric value, a step for the grid of h for continuous variable.
#' @param p1_objmethod a character string, the method using for the tuning.
#' @param p2_objmethod a character string, the method using for the tuning.
#' @param p1_lat_X_min a positive integer, minimum size of the latent space of the common variables in phase 1.
#' @param p1_lat_X_max a positive integer, maximum size of the latent space of the common variables in phase 1.
#' @param p1_lat_X_step a positive integer, step size of the latent space of the common variables in phase 1.
#' @param p2_lat_X_min a positive integer, minimum size of the latent space of the common variables in phase 2.
#' @param p2_lat_X_max a positive integer, maximum size of the latent space of the common variables in phase 2.
#' @param p2_lat_X_step a positive integer, step size of the latent space of the common variables in phase 2.
#' @param p1_lat_Y_min a positive integer, minimum size of the latent space of the non-common variables in phase 1.
#' @param p1_lat_Y_max a positive integer, maximum size of the latent space of the non-common variables in phase 1.
#' @param p1_lat_Y_step a positive integer, step size of the latent space of the non-common variables in phase 1.
#' @param p2_lat_Y_min a positive integer, minimum size of the latent space of the non-common variables in phase 2.
#' @param p2_lat_Y_max a positive integer, maximum size of the latent space of the non-common variables in phase 2.
#' @param p2_lat_Y_step a positive integer, step size of the latent space of the non-common variables in phase 2.
#' @param p1_n_combs a positive integer, number of random combinations of tuning parameters generated in phase 1 (used only when \code{tuning_type = "random"}).
#' @param p2_n_combs a positive integer, number of random combinations of tuning parameters generated in phase 2 (used only when \code{tuning_type = "random"}).
#' @param p1_man_params a data.frame, a manual way to specify the grid of tuning parameters in phase 1.
#' @param p2_man_params a data.frame, a manual way to specify the grid of tuning parameters in phase 2.
#' @param exactMatch a logical, indicate if an exact matching has to be performed or not.
#' @return A list with all the options required by the statistical matching algorithm based on PCA-CCA.
#' @details DETAILS
#'
#'
#' @export

statMatch.PCA.options <- function(print.details = TRUE,
                                  print.details_tuning = FALSE,
                                  scaling = c("z-score", "min-max", "no"),
                                  tuning_type = c("grid", "random"),
                                  type_predict = c("loop", "matrix"),
                                  exactMatch = FALSE,
                                  names.CV_Cat_prio = NULL,
                                  n_fold = 5L,
                                  # Arguments for CCA
                                  #d = 2L,
                                  d = NULL,
                                  p1_dmax = 2L,
                                  p1_dmin = 2L,
                                  p1_dstep = 2L,
                                  p2_dmax = 2L,
                                  p2_dmin = 2L,
                                  p2_dstep = 2L,
                                  rot = FALSE,
                                  p1_kernel_predict = c("alea", "epan", "gauss"),
                                  p1_kernel_predict_tuning = NULL,
                                  p2_kernel_predict = c("gauss", "epan", "alea"),
                                  # Arguments for PCA
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
                                  # Argument for the tuning
                                  p1_objmethod = c("wRMSE", "wsRMSE", "wMCR", "wCE"),
                                  p2_objmethod = c("wsRMSE", "wRMSE", "wCE"),
                                  p1_h = NULL,
                                  p1_hmax = 1,
                                  p1_hmin = 0.01,
                                  p1_hstep = 0.1,
                                  p2_h = NULL,
                                  p2_hmax = 1,
                                  p2_hmin = 0.01,
                                  p2_hstep = 0.1,
                                  p1_n_combs = 10L,
                                  p2_n_combs = 10L,
                                  p1_man_params = NULL,
                                  p2_man_params = NULL) {
  scaling <- match.arg(scaling)
  tuning_type <- match.arg(tuning_type)
  type_predict <- match.arg(type_predict)

  # Check input variables
  p1_kernel_predict <- match.arg(p1_kernel_predict)
  if (is.null(p1_kernel_predict_tuning)) p1_kernel_predict_tuning <- p1_kernel_predict
  p2_kernel_predict <- match.arg(p2_kernel_predict)
  p1_objmethod <- match.arg(p1_objmethod)
  p2_objmethod <- match.arg(p2_objmethod)
  if (length(n_fold) != 1 | !is.integer(n_fold) | (n_fold <= 0))
    stop("PCA.options: argument 'n_fold' should be a positive integer.")

  if (length(p1_dmax) != 1 | !is.integer(p1_dmax) |  (p1_dmax <= 0))
    stop("KCCA.options: argument 'p1_dmax' should be a positive integer")

  if (length(p1_dmin) != 1 | !is.integer(p1_dmin) |  (p1_dmin <= 0))
    stop("KCCA.options: argument 'p1_dmax' should be a positive integer")

  if ((p1_dmax >= 3) & (p1_kernel_predict == "epan"))
    stop("KCCA.options: 'p1_dmax' bigger than two is not implement with epan")

  if ((p2_dmax >= 3) & (p2_kernel_predict == "epan"))
    stop("KCCA.options: 'd' bigger than two is not implement with epan")

  if (length(exactMatch) != 1 | !is.logical(exactMatch) | any(is.na(exactMatch)))
    stop("superPCA.options: argument 'exactMatch' should be a logical.")


  if (length(rot) != 1 | !is.logical(rot) | any(is.na(rot)))
    stop("superPCA.options: argument 'rot' should be a logical.")

  if ((p1_kernel_predict == "epan") & (type_predict == "loop"))
    stop("PCA.options: epan kernel is only possible with matrix")

  if ((p2_kernel_predict == "epan") & (type_predict == "loop"))
    stop("PCA.options: epan kernel is only possible with matrix")

  if ((p1_kernel_predict == "gauss") & (type_predict == "loop"))
    stop("PCA.options: gauss kernel in phase 1 is only possible with matrix")

  if ((p2_kernel_predict == "alea") & (type_predict == "loop"))
    stop("PCA.options: alea kernel in phase 2 is only possible with matrix")

  if (is.null(d)) {
    p1_d <- seq(p1_dmin, p1_dmax, by = p1_dstep)
    p2_d <- seq(p2_dmin, p2_dmax, by = p2_dstep)
  } else {
    p1_d <- d
    p2_d <- d
  }

  if (is.null(p1_h)) p1_h <- seq(p1_hmin, p1_hmax, by = p1_hstep)
  if (is.null(p2_h)) p2_h <- seq(p2_hmin, p2_hmax, by = p2_hstep)

  if (is.null(p1_man_params)) {
    if (length(p1_n_combs) != 1 | !is.integer(p1_n_combs) | (p1_n_combs <= 0))
      stop("PCA.options: argument 'p1_n_combs' should be a positive integer.")

    if (length(p1_lat_X_min) != 1 | !is.integer(p1_lat_X_min) | (p1_lat_X_min <= 0))
      stop("PCA.options: argument 'p1_lat_X_min' should be a positive integer.")

    if (length(p1_lat_X_max) != 1 | !is.integer(p1_lat_X_max) | (p1_lat_X_max <= 0) | (p1_lat_X_max < p1_lat_X_min))
      stop("PCA.options: argument 'p1_lat_X_max' should be a positive integer greater or equal to 'p1_lat_X_min'.")

    if (length(p1_lat_X_step) != 1 | !is.integer(p1_lat_X_step) | (p1_lat_X_step <= 0))
      stop("PCA.options: argument 'p1_lat_X_step' should be a positive integer.")

    if (length(p1_lat_Y_min) != 1 | !is.integer(p1_lat_Y_min) | (p1_lat_Y_min <= 0))
      stop("PCA.options: argument 'p1_lat_Y_min' should be a positive integer.")

    if (length(p1_lat_Y_max) != 1 | !is.integer(p1_lat_Y_max) | (p1_lat_Y_max <= 0) | (p1_lat_Y_max < p1_lat_Y_min))
      stop("PCA.options: argument 'p1_lat_Y_max' should be a positive integer greater or equal to 'p1_lat_Y_min'.")

    if (length(p1_lat_Y_step) != 1 | !is.integer(p1_lat_Y_step) | (p1_lat_Y_step <= 0))
      stop("PCA.options: argument 'p1_lat_Y_step' should be a positive integer.")

    if (length(p1_h) == 0 | !is.numeric(p1_h) | any(p1_h <= 0))
      stop("PCA.options: argument 'p1_h' should be a numeric vector containing double(s) greater than zero.")

    if (length(p1_d) == 0 | !is.numeric(p1_d) | any(p1_d <= 0))
      stop("KCCA.options: argument 'p1_d' should be a numeric vector containing double(s) greater than zero.")

    if (p1_lat_X_min < max(p1_d) | p1_lat_Y_min < max(p1_d))
      stop("PCA.options: arguments 'p1_lat_X_min', 'p1_lat_Y_min' should be greater or equal to 'max(p1_d)' (the number of latent variables in CCA).")

    hpar1 <- init_PCA(
      p1_h, p1_lat_X_min, p1_lat_X_max, p1_lat_X_step,
      p1_lat_Y_min, p1_lat_Y_max, p1_lat_Y_step, p1_d, tuning_type, p1_n_combs
    )
  } else {
    hpar1 <- check_man_params_PCA(p1_man_params, "1")
    p1_h <- unique(hpar1$h)
  }

  if (is.null(p2_man_params)) {
    if (length(p2_n_combs) != 1 | !is.integer(p2_n_combs) | (p2_n_combs <= 0))
      stop("PCA.options: argument 'p2_n_combs' should be a positive integer.")

    if (length(p2_lat_X_min) != 1 | !is.integer(p2_lat_X_min) | (p2_lat_X_min <= 0))
      stop("PCA.options: argument 'p2_lat_X_min' should be a positive integer.")

    if (length(p2_lat_X_max) != 1 | !is.integer(p2_lat_X_max) | (p2_lat_X_max <= 0) | (p2_lat_X_max < p2_lat_X_min))
      stop("PCA.options: argument 'p2_lat_X_max' should be a positive integer greater or equal to 'p2_lat_X_min'.")

    if (length(p2_lat_X_step) != 1 | !is.integer(p2_lat_X_step) | (p2_lat_X_step <= 0))
      stop("PCA.options: argument 'p2_lat_X_step' should be a positive integer.")

    if (length(p2_lat_Y_min) != 1 | !is.integer(p2_lat_Y_min) | (p2_lat_Y_min <= 0))
      stop("PCA.options: argument 'p2_lat_Y_min' should be a positive integer.")

    if (length(p2_lat_Y_max) != 1 | !is.integer(p2_lat_Y_max) | (p2_lat_Y_max <= 0) | (p2_lat_Y_max < p2_lat_Y_min))
      stop("PCA.options: argument 'p2_lat_Y_max' should be a positive integer greater or equal to 'p2_lat_Y_min'.")

    if (length(p2_lat_Y_step) != 1 | !is.integer(p2_lat_Y_step) | (p2_lat_Y_step <= 0))
      stop("PCA.options: argument 'p2_lat_Y_step' should be a positive integer.")

    if (length(p2_h) == 0 | !is.numeric(p2_h) | any(p2_h <= 0))
      stop("PCA.options: argument 'p2_h' should be a numeric vector containing double(s) greater than zero.")

    if (length(p2_d) == 0 | !is.numeric(p2_d) | any(p2_d <= 0))
      stop("KCCA.options: argument 'p1_d' should be a numeric vector containing double(s) greater than zero.")

    if (p2_lat_X_min < max(p2_d) | p2_lat_Y_min < max(p2_d))
      stop("PCA.options: arguments 'p2_lat_X_min', 'p2_lat_Y_min' should be greater or equal to 'max(p2_d)' (the number of latent variables in CCA).")

    hpar2 <- init_PCA(
      p2_h, p2_lat_X_min, p2_lat_X_max, p2_lat_X_step,
      p2_lat_Y_min, p2_lat_Y_max, p2_lat_Y_step, p2_d, tuning_type, p2_n_combs
    )
  } else {
    hpar2 <- check_man_params_PCA(p2_man_params, "2")
    p2_h <- unique(hpar2$h)
  }

  if (tuning_type == "grid") {
    p1_n_combs <- nrow(hpar1) / length(p1_h)
    p2_n_combs <- nrow(hpar2) / length(p2_h)
  }

  P1 <- list(
    n_combs = p1_n_combs,
    hpar = hpar1,
    n_h = length(p1_h),
    kernel_predict = p1_kernel_predict,
    kernel_predict_tuning = p1_kernel_predict_tuning,
    objmethod = p1_objmethod
  )

  P2 <- list(
    n_combs = p2_n_combs,
    hpar = hpar2,
    n_h = length(p2_h),
    kernel_predict = p2_kernel_predict,
    kernel_predict_tuning = p2_kernel_predict,
    objmethod = p2_objmethod
  )

  opts <- list(
    print.details = print.details,
    print.details_tuning = print.details_tuning,
    scaling = scaling,
    n_fold = n_fold,
    exactMatch = exactMatch,
    names.CV_Cat_prio = names.CV_Cat_prio,
    n_fold = n_fold,

    # Arguments for CCA
    #d = d,
    rot = rot,
    type_predict = type_predict,

    # Argument for the tuning
    P1 = P1,
    P2 = P2
  )
  class(opts) <- "PCA.options"
  return(opts)
}



init_PCA <- function(h, lat_X_min, lat_X_max, lat_X_step, lat_Y_min, lat_Y_max, lat_Y_step, d, tuning_type, n_combs) {
  if (tuning_type == "grid") {
    out <- merge(as.data.frame(h), as.data.frame(seq(from = lat_X_min, to = lat_X_max, by = lat_X_step)))
    out <- merge(out, as.data.frame(seq(from = lat_Y_min, to = lat_Y_max, by = lat_Y_step)))
    out <- merge(out, as.data.frame(d))
  } else if (tuning_type == "random") {
    out <- data.frame(
      h = sample(h, size = n_combs, replace = TRUE),
      lat_X = sample(seq(from = lat_X_min, to = lat_X_max, by = lat_X_step), size = n_combs, replace = TRUE),
      lat_Y = sample(seq(from = lat_Y_min, to = lat_Y_max, by = lat_Y_step), size = n_combs, replace = TRUE),
      d = sample(d, size = n_combs, replace = TRUE)
    )
  }

  names(out) <- c("h", "lat_X", "lat_Y", "d")
  return(out)
}


check_man_params_PCA <- function(params, phase = c("1", "2")) {
  phase <- match.arg(phase)
  if (!is.data.frame(params)) {
    stop("PCA.options: argument 'man_params_", phase, "' should be a data.frame object.")
  }
  # Check if there is any missing column
  req_names <- c("h", "lat_X", "lat_Y", "d")
  check <- !req_names %in% colnames(params)
  missing <- paste0("'", req_names[check], "'", collapse = ", ")
  if (any(check)) {
    stop("PCA.options: column(s) ", missing, " missing in 'man_params_", phase, "' data.frame.")
  }
  # Check numeric format
  check_names <- c("h", "lat_X", "lat_Y",, "d")
  for (i in 1:length(check_names)) {
    if (!is.numeric(params[, check_names[i]]) | any(params[, check_names[i]] <= 0)) {
      stop(
        "PCA.options: column '", check_names[i], "' in 'man_params_", phase,
        "' data.frame should be a numeric vector with values greater than zero."
      )
    }
  }
  if( nrow(params) == length(unique(params$h))*length(unique(params$hx))*length(unique(params$hy))*length(unique(params$g))*length(unique(params$d)))
  { return(params)
  } else {
    if( nrow(params) %%length(unique(params$h)) == 0)
    {
      return(params)
    } else {
      stop("number h must be the same for each triplet of hx, hy g and d.")
    }
  }
}
