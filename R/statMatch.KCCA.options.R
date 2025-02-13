#' @title Statistical Matching: KCCA options
#' @description Function to create the options for KCCA.
#'
#' @param print.details Logical Value, print details or not (in the tuning part only print the value of the objective function for each parameters combination).
#' @param print.details_tuning Logical Value, if true print all the details in the tuning part.
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
#' @param p1_h a numeric vector, vector of value for h, the bandwidth of the Matching kernel for categorical variable.
#' @param p1_hmax a numeric value, maximal value for h, the bandwidth of the Matching kernel for categorical variable.
#' @param p1_hmin a numeric value, minimal value for h, the bandwidth of the Matching kernel for categorical variable.
#' @param p1_hstep a numeric value, step for the grid of h for categorical variable.
#' @param p2_h a numeric vector, vector of value for h, the bandwidth of the Matching kernel for continuous variable.
#' @param p2_hmax a numeric value, maximal value for h, the bandwidth of the Matching kernel for continuous variable.
#' @param p2_hmin a numeric value, minimal value for h, the bandwidth of the Matching kernel for continuous variable.
#' @param p2_hstep a numeric value, step for the grid of h for continuous variable.
#' @param p1_hx a numeric vector, vector of value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p1_hxmax  a numeric value, maximal value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p1_hxmin  a numeric value, minimal value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p1_hxstep  a numeric value,  step for the grid of h for categorical variable.
#' @param p2_hx  a numeric vector, vector of value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p2_hxmax  a numeric value, maximal value for h, the bandwidth of the KCCA kernel for continuous variable.
#' @param p2_hxmin  a numeric value, minimal value for h, the bandwidth of the KCCA kernel for continuous variable.
#' @param p2_hxstep  a numeric value, step for the grid of h for continuous variable.
#' @param p1_hy a numeric vector, vector of value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p1_hymax  a numeric value, maximal value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p1_hymin  a numeric value, minimal value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p1_hystep  a numeric value, step for the grid of h for categorical variable.
#' @param p2_hy a numeric vector, vector of value for h, the bandwidth of the KCCA kernel for continuous variable.
#' @param p2_hymax  a numeric value, maximal value for h, the bandwidth of the KCCA kernel for continuous variable.
#' @param p2_hymin  a numeric value, minimal value for h, the bandwidth of the KCCA kernel for continuous variable.
#' @param p2_hystep  a numeric value, step for the grid of h for continuous variable.
#' @param p1_g a numeric vector, vector of value for g, the regularization parameter for categorical variable.
#' @param p1_gmax  a numeric value, maximal value for g, the regularization parameter for categorical variable.
#' @param p1_gmin  a numeric value, minimal value for g, the regularization parameter for categorical variable.
#' @param p1_gstep  a numeric value, mtep for the grid of g for categorical variable.
#' @param p2_g a numeric vector, vector of value for g, the regularization parameter for continuous variable.
#' @param p2_gmax  a numeric value, maximal value for g, the regularization parameter for continuous variable.
#' @param p2_gmin  a numeric value, minimal value for g, the regularization parameter for continuous variable.
#' @param p2_gstep  a numeric value, step for the grid of g for continuous variable.
#' @param p1_objmethod a character string, Method using for the tuning.
#' @param p2_objmethod a character string, Method using for the tuning.
#' @param p1_n_combs a positive integer, number of random combinations of tuning parameters generated in phase 1 (used only when \code{tuning_type = "random"}).
#' @param p2_n_combs a positive integer, number of random combinations of tuning parameters generated in phase 2 (used only when \code{tuning_type = "random"}).
#' @param p1_man_params a data.frame, a manual way to specify the grid of tuning parameters in phase 1.
#' @param p2_man_params a data.frame, a manual way to specify the grid of tuning parameters in phase 2.
#' @param nc a positive integer, number of core to use during the tuning phase
#' @param par a logical, if TRUE doing parallelisation during the tuning phase
#' @param exactMatch a logical, indicate if an exact matching has to be performed or not.
#' @return A list with all the options required by the statistical matching algorithm based on KCCA.
#' @details DETAILS
#' @export

statMatch.KCCA.options <- function(print.details = TRUE,
                                   print.details_tuning = FALSE,
                                   scaling = c("z-score", "min-max", "no"),
                                   tuning_type = c("two h","one h","random"),
                                   type_predict =c("loop","matrix"),
                                   exactMatch = FALSE,
                                   names.CV_Cat_prio = NULL,
                                   n_fold = 5L,
                                   # Arguments for parallelisation
                                   par = FALSE,
                                   par_fold = TRUE,
                                   nc = 2L,
                                   # Arguments for CCA or KCCA
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
                                   # Argument for the tuning
                                   p1_objmethod = c("wRMSE","wsRMSE", "wMCR", "wCE"),
                                   p2_objmethod = c("wsRMSE", "wRMSE", "wCE"),
                                   p1_h = NULL,
                                   p1_hmax = 1,
                                   p1_hmin = 0.01,
                                   p1_hstep = 0.1,
                                   p1_hx = NULL,
                                   p1_hxmax = 1,
                                   p1_hxmin = 0.01,
                                   p1_hxstep = 0.25,
                                   p1_hy = NULL,
                                   p1_hymax = NULL,
                                   p1_hymin = NULL,
                                   p1_hystep  = NULL,
                                   p1_g = 0.00002,
                                   p1_gmax = 0.0001,
                                   p1_gmin = 0.00001,
                                   p1_gstep = 0.00001,
                                   p2_h = NULL,
                                   p2_hmax = 1,
                                   p2_hmin = 0.01,
                                   p2_hstep = 0.1,
                                   p2_hx = NULL,
                                   p2_hxmax = 1,
                                   p2_hxmin = 0.01,
                                   p2_hxstep = 0.25,
                                   p2_hy = NULL,
                                   p2_hymax = NULL,
                                   p2_hymin = NULL,
                                   p2_hystep = NULL,
                                   p2_g = 0.00002,
                                   p2_gmax = 0.0001,
                                   p2_gmin = 0.00001,
                                   p2_gstep = 0.00001,
                                   p1_n_combs = 10L,
                                   p2_n_combs = 10L,
                                   p1_man_params = NULL,
                                   p2_man_params = NULL) {

  scaling <- match.arg(scaling)
  tuning_type <- match.arg(tuning_type)
  type_predict <- match.arg(type_predict)

  # Check input variables
  p1_kernel_predict <- match.arg(p1_kernel_predict)
  if (is.null(p1_kernel_predict_tuning))  p1_kernel_predict_tuning <- p1_kernel_predict
  p2_kernel_predict <- match.arg(p2_kernel_predict)
  p1_objmethod <- match.arg(p1_objmethod)
  p2_objmethod <- match.arg(p2_objmethod)
  if (length(n_fold) != 1 | !is.integer(n_fold) | (n_fold <= 0))
    stop("KCCA.options: argument 'n_fold' should be a positive integer.")

  if (length(nc) != 1 | !is.integer(nc) | (nc <= 0))
    stop("KCCA.options: argument 'nc' should be a positive integer.")

  if (length(p1_dmax) != 1 | !is.integer(p1_dmax) |  (p1_dmax <= 0))
    stop("KCCA.options: argument 'p1_dmax' should be a positive integer")

  if (length(p1_dmin) != 1 | !is.integer(p1_dmin) |  (p1_dmin <= 0))
    stop("KCCA.options: argument 'p1_dmax' should be a positive integer")

 if ((p1_dmax >= 3) & (p1_kernel_predict == "epan"))
   stop("KCCA.options: 'p1_dmax' bigger than two is not implement with epan")

  if ((p2_dmax >= 3) & (p2_kernel_predict == "epan"))
    stop("KCCA.options: 'd' bigger than two is not implement with epan")

  if (length(exactMatch) != 1 | !is.logical(exactMatch) | any(is.na(exactMatch)))
    stop("superKCCA.options: argument 'exactMatch' should be a logical.")

  if (length(rot) != 1 | !is.logical(rot) | any(is.na(rot)))
    stop("superKCCA.options: argument 'rot' should be a logical.")

  if (length(par) != 1 | !is.logical(par) | any(is.na(par)))
    stop("superKCCA.options: argument 'par' should be a logical.")

  if (length(par_fold) != 1 | !is.logical(par_fold) | any(is.na(par_fold)))
    stop("superKCCA.options: argument 'par_fold' should be a logical.")

  if ((p1_kernel_predict == "epan") & (type_predict == "loop"))
    stop("KCCA.options: epan kernel is only possible with matrix")

  if ((p2_kernel_predict == "epan") & (type_predict == "loop"))
    stop("KCCA.options: epan kernel is only possible with matrix")

  if ((p1_kernel_predict == "gauss") & (type_predict == "loop"))
    stop("KCCA.options: gauss kernel in phase 1 is only possible with matrix")

  if ((p2_kernel_predict == "alea") & (type_predict == "loop"))
    stop("KCCA.options: alea kernel in phase 2 is only possible with matrix")

  if (is.null(p1_hymax)) p1_hymax <- p1_hxmax
  if (is.null(p1_hymin)) p1_hymin <- p1_hxmin
  if (is.null(p1_hystep)) p1_hystep <- p1_hxstep
  if (is.null(p2_hymax)) p2_hymax <- p2_hxmax
  if (is.null(p2_hymin)) p2_hymin <- p2_hxmin
  if (is.null(p2_hystep)) p2_hystep <- p2_hxstep

  if (is.null(d)) {
    p1_d <- seq(p1_dmin, p1_dmax, by = p1_dstep)
    p2_d <- seq(p2_dmin, p2_dmax, by = p2_dstep)
  } else {
    p1_d <- d
    p2_d <- d
  }
  if (is.null(p1_h)) p1_h <- seq(p1_hmin, p1_hmax, by = p1_hstep)
  if (is.null(p1_hx)) p1_hx <- seq(p1_hxmin, p1_hxmax, by = p1_hxstep)
  if (is.null(p1_hy)) p1_hy <- seq(p1_hymin, p1_hymax, by = p1_hystep)
  if (is.null(p1_g)) p1_g <- seq(p1_gmin, p1_gmax, by = p1_gstep)

  if (is.null(p2_h)) p2_h <- seq(p2_hmin, p2_hmax, by = p2_hstep)
  if (is.null(p2_hx)) p2_hx <- seq(p2_hxmin, p2_hxmax, by = p2_hxstep)
  if (is.null(p2_hy)) p2_hy <- seq(p2_hymin, p2_hymax, by = p2_hystep)
  if (is.null(p2_g)) p2_g <- seq(p2_gmin, p2_gmax, by = p2_gstep)


  if (tuning_type == "two h") {
    p1_n_combs <- length(p1_hx)*length(p1_hy)*length(p1_g)*length(p1_d)
    p2_n_combs <- length(p2_hx)*length(p2_hy)*length(p2_g)*length(p1_d)
  } else if (tuning_type == "one h") {
    p1_n_combs <- length(p1_hx)*length(p1_g)*length(p1_d)
    p2_n_combs <- length(p2_hx)*length(p2_g)*length(p1_d)
  }

  if (is.null(p1_man_params)) {
    if (length(p1_n_combs) != 1 | !is.integer(p1_n_combs) | (p1_n_combs <= 0))
      stop("KCCA.options: argument 'p1_n_combs' should be a positive integer.")
    if (length(p1_h) == 0 | !is.numeric(p1_h) | any(p1_h <= 0))
      stop("KCCA.options: argument 'p1_h' should be a numeric vector containing double(s) greater than zero.")
    if (length(p1_hx) == 0 | !is.numeric(p1_hx) | any(p1_h <= 0))
      stop("KCCA.options: argument 'p1_hx' should be a numeric vector containing double(s) greater than zero.")
    if (length(p1_hy) == 0 | !is.numeric(p1_hy) | any(p1_hy <= 0))
      stop("KCCA.options: argument 'p1_hy' should be a numeric vector containing double(s) greater than zero.")
    if (length(p1_g) == 0 | !is.numeric(p1_g) | any(p1_g <= 0))
      stop("KCCA.options: argument 'p1_g' should be a numeric vector containing double(s) greater than zero.")
    if (length(p1_d) == 0 | !is.numeric(p1_d) | any(p1_d <= 0))
      stop("KCCA.options: argument 'p1_d' should be a numeric vector containing double(s) greater than zero.")
     hpar1 <- init_KCCA(p1_h,p1_hx,p1_hy,p1_g,p1_d,tuning_type)
  } else {
    hpar1 <- check_man_params_KCCA(p1_man_params,"1")
    p1_h <- unique(hpar1$h)
    p1_n_combs <- nrow(hpar1)/length(p1_h)
  }

  if (is.null(p2_man_params)) {
    if (length(p2_n_combs) != 1 | !is.integer(p2_n_combs) | (p2_n_combs <= 0))
      stop("KCCA.options: argument 'p2_n_combs' should be a positive integer.")
    if (length(p2_h) == 0 | !is.numeric(p2_h) | any(p2_h <= 0))
      stop("KCCA.options: argument 'p2_h' should be a numeric vector containing double(s) greater than zero.")
    if (length(p2_hx) == 0 | !is.numeric(p2_hx) | any(p2_h <= 0))
      stop("KCCA.options: argument 'p2_hx' should be a numeric vector containing double(s) greater than zero.")
    if (length(p2_hy) == 0 | !is.numeric(p2_hy) | any(p2_hy <= 0))
      stop("KCCA.options: argument 'p2_hy' should be a numeric vector containing double(s) greater than zero.")
    if (length(p2_g) == 0 | !is.numeric(p2_g) | any(p2_g <= 0))
      stop("KCCA.options: argument 'p2_g' should be a numeric vector containing double(s) greater than zero.")
    if (length(p2_d) == 0 | !is.numeric(p2_d) | any(p2_d <= 0))
      stop("KCCA.options: argument 'p2_d' should be a numeric vector containing double(s) greater than zero.")


   hpar2 <- init_KCCA(p2_h,p2_hx,p2_hy,p2_g,p2_d,tuning_type)
  } else {
    hpar2 <- check_man_params_KCCA(p2_man_params,"2")
    p2_h <- unique(hpar2$h)
    p2_n_combs <- nrow(hpar2)/length(p2_h)
  }

  P1 <- list(
    #d = p1_d,
    n_combs = p1_n_combs ,
    hpar = hpar1,
    n_h = length(p1_h),
    kernel_predict = p1_kernel_predict,
    kernel_predict_tuning = p1_kernel_predict_tuning,
    objmethod = p1_objmethod
  )

  P2 <- list(
    #d = p2_d,
    n_combs = p2_n_combs,
    hpar = hpar2,
    n_h = length(p2_h),
    kernel_predict=p2_kernel_predict,
    kernel_predict_tuning=p2_kernel_predict,
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

    # Arguments for CCA or KCCA
    rot = rot,
    type_predict = type_predict,

    # Argument for the tuning
    P1 = P1,
    P2 = P2,
    par = par,
    par_fold = par_fold,
    nc = nc
  )
  class(opts) <- "KCCA.options"
  return(opts)
}



init_KCCA <- function(h,hx,hy,g,d,tuning_type,n_combs) {

  if (tuning_type=="two h") {
  out <- merge(as.data.frame(h),as.data.frame(hx))
  out <- merge(out,as.data.frame(hy))
  out <- merge(out,as.data.frame(g))
  out <- merge(out,as.data.frame(d))
  } else if (tuning_type=="one h") {
    out <- data.frame(hx=hx,hy=hx) #Same h for x and y
    out <- merge(as.data.frame(h),out)
    out <- merge(out,as.data.frame(g))
    out <- merge(out,as.data.frame(d))
  } else if (tuning_type=="random") {
    out <- data.frame(
      hx = sample(hx, size = n_combs, replace = TRUE),
      hy = sample(hy, size = n_combs, replace = TRUE),
      g = sample(g, size = n_combs, replace = TRUE))

    out <- merge(as.data.frame(h), out)
    out <- merge(out,as.data.frame(d))
  }

  names(out) <- c("h","hx","hy","g", "d")
  return(out)
}


check_man_params_KCCA <- function(params, phase = c("1", "2")) {
  phase <- match.arg(phase)
  if (!is.data.frame(params))
    stop("KCCA.options: argument 'man_params_", phase, "' should be a data.frame object.")
  # Check if there is any missing column
  req_names <- c("h", "hx","hy","g", "d")
  check <- !req_names %in% colnames(params)
  missing <- paste0("'", req_names[check],"'", collapse = ", ")
  if (any(check))
    stop("KCCA.options: column(s) ", missing, " missing in 'man_params_", phase, "' data.frame.")
  # Check numeric format
  check_names <-  c("h", "hx","hy","g", "d")
  for (i in 1:length(check_names)) {
    if (!is.numeric(params[,check_names[i]]) | any(params[,check_names[i]] <= 0))
      stop("KCCA.options: column '", check_names[i], "' in 'man_params_", phase,
           "' data.frame should be a numeric vector with values greater than zero.")
  }
   if( nrow(params) == length(unique(params$h))*length(unique(params$hx))*length(unique(params$hy))*length(unique(params$g))*length(unique(params$d)))
   { return(params)
   } else {
     if( nrow(params) %%length(unique(params$h)) == 0)
         {
           return(params)
         } else {
           stop("number h must be the same for each triplet of hx, hy and g")
         }
     }

}
