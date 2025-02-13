#' @title Statistical Matching: KCCA options
#' @description Function to create the options for KCCA.
#'
#' @param print.details Logical Value, print details or not (in the tuning part only print the value of the objective function for each parameters combination).
#' @param print.details_tuning Logical Value, if true print all the details in the tuning part.
#' @param tuning_type a character string, the method used to find the optimal tuning parameters: \code{"random"}, \code{"two h"} or \code{"one h"}
#' @param scaling a character string, the name of the method used to scale the variables in the receiver and donor data sets: \code{"z-score"}, \code{"min-max"} or \code{"no"} (no scaling).
#' @param d a positive integer, number of latent variable used in CCA.
#' @param p0_dmax a positive integer, number of latent variable used in CCA.
#' @param p0_dmin a positive integer, number of latent variable used in CCA.
#' @param p0_dstep a positive integer, number of latent variable used in CCA.
#' @param rot a logical, if TRUE the bandwidth is h multiply by the variance.
#' @param type_predict a character string, Type of prediction \code{"matrix"} or \code{"loop"}.
#' @param n_fold a positive integer, number of folds used for the cross-validation of the tuning parameters.
#' @param names.CV_Cat_prio a character vector, Name of common categorical variable that have to match.
#' @param p0_kernel_predict a character string, Type of kernel use for prediction for continuous variable.
#' @param p0_kernel_predict_cat a character string, Type of kernel use for prediction for categorical variable.
#' @param p0_h a numeric vector, vector of value for h, the bandwidth of the Matching kernel for categorical variable.
#' @param p0_hmax a numeric value, maximal value for h, the bandwidth of the Matching kernel for categorical variable.
#' @param p0_hmin a numeric value, minimal value for h, the bandwidth of the Matching kernel for categorical variable.
#' @param p0_hstep a numeric value, step for the grid of h for categorical variable.
#' @param p0_hx a numeric vector, vector of value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p0_hxmax  a numeric value, maximal value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p0_hxmin  a numeric value, minimal value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p0_hxstep  a numeric value,  step for the grid of h for categorical variable.
#' @param p0_hy a numeric vector, vector of value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p0_hymax  a numeric value, maximal value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p0_hymin  a numeric value, minimal value for h, the bandwidth of the KCCA kernel for categorical variable.
#' @param p0_hystep  a numeric value, step for the grid of h for categorical variable.
#' @param p0_g a numeric vector, vector of value for g, the regularization parameter for categorical variable.
#' @param p0_gmax  a numeric value, maximal value for g, the regularization parameter for categorical variable.
#' @param p0_gmin  a numeric value, minimal value for g, the regularization parameter for categorical variable.
#' @param p0_gstep  a numeric value, mtep for the grid of g for categorical variable.
#' @param p0_objmethod a character string, Method using for the tuning.
#' @param p0_n_combs a positive integer, number of random combinations of tuning parameters generated in phase 1 (used only when \code{tuning_type = "random"}).
#' @param p0_man_params a data.frame, a manual way to specify the grid of tuning parameters in phase 1.
#' @param nc a positive integer, number of core to use during the tuning phase
#' @param par a logical, if TRUE doing parallelisation during the tuning phase
#' @param exactMatch a logical, indicate if an exact matching has to be performed or not.
#' @return A list with all the options required by the statistical matching algorithm based on KCCA.
#' @details DETAILS
#' @export

statMatch.KCCA.one.phase.options <- function(print.details = TRUE,
                                   print.details_tuning = FALSE,
                                   scaling = c("z-score", "min-max", "no"),
                                   tuning_type = c("two h","one h","random"),
                                   tuning_type_bandwith = c("one h","two h"),
                                   type_predict =c("loop","matrix"),
                                   exactMatch = FALSE,
                                   names.CV_Cat_prio = NULL,
                                   n_fold = 5L,
                                   # Arguments for parallelisation
                                   par = FALSE,
                                   par_fold = FALSE,
                                   nc = 2L,
                                   # Arguments for CCA or KCCA
                                   d = NULL,
                                   p0_dmax = 2L,
                                   p0_dmin = 2L,
                                   p0_dstep = 2L,
                                   rot = FALSE,
                                   p0_kernel_predict = c("gauss", "alea", "epan"),
                                   p0_kernel_predict_cat = c("alea", "epan", "gauss"),
                                   # Argument for the tuning
                                   p0_objmethod = c("wsRMSE", "wRMSE", "wMCR", "wCE"),
                                   p0_h = NULL,
                                   p0_hmax = 1,
                                   p0_hmin = 0.01,
                                   p0_hstep = 0.1,
                                   p0_h_cat = NULL,
                                   p0_hmax_cat = NULL,
                                   p0_hmin_cat = NULL,
                                   p0_hstep_cat = NULL,
                                   p0_hx = NULL,
                                   p0_hxmax = 1,
                                   p0_hxmin = 0.01,
                                   p0_hxstep = 0.25,
                                   p0_hy = NULL,
                                   p0_hymax = NULL,
                                   p0_hymin = NULL,
                                   p0_hystep  = NULL,
                                   p0_g = NULL,
                                   p0_gmax = 0.5,
                                   p0_gmin = 0.01,
                                   p0_gstep = 0.1,
                                   p0_n_combs = 10L,
                                   p0_man_params = NULL) {

  scaling <- match.arg(scaling)
  tuning_type <- match.arg(tuning_type)
  tuning_type_bandwith <- match.arg(tuning_type_bandwith)
  type_predict <- match.arg(type_predict)

  # Check input variables
  p0_kernel_predict <- match.arg(p0_kernel_predict)
  p0_kernel_predict_cat <- match.arg(p0_kernel_predict_cat)
  p0_objmethod <- match.arg(p0_objmethod)
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

  if ((p0_kernel_predict == "epan") & (type_predict == "loop"))
    stop("KCCA.options: epan kernel is only possible with matrix")


  if (is.null(p0_hymax)) p0_hymax <- p0_hxmax
  if (is.null(p0_hymin)) p0_hymin <- p0_hxmin
  if (is.null(p0_hystep)) p0_hystep <- p0_hxstep


  if (is.null(p0_hmax_cat)) p0_hmax_cat <- p0_hmax
  if (is.null(p0_hmin_cat)) p0_hmin_cat <- p0_hmin
  if (is.null(p0_hstep_cat)) p0_hstep_cat <- p0_hstep

  if (is.null(p0_h)) p0_h <- seq(p0_hmin, p0_hmax, by = p0_hstep)
  if (is.null(p0_hx)) p0_hx <- seq(p0_hxmin, p0_hxmax, by = p0_hxstep)
  if (is.null(p0_hy)) p0_hy <- seq(p0_hymin, p0_hymax, by = p0_hystep)
  if (is.null(p0_g)) p0_g <- seq(p0_gmin, p0_gmax, by = p0_gstep)

  if (is.null(d)) {
    p0_d <- seq(p0_dmin, p0_dmax, by = p0_dstep)
  } else {
    p0_d <- d
  }

  if (is.null(p0_h_cat)) {
    if (tuning_type_bandwith == "two h") {
      p0_h_cat <- seq(p0_hmin_cat, p0_hmax_cat, by = p0_hstep_cat)
    } else if (tuning_type_bandwith == "one h") {
      p0_h_cat <- p0_h
    }
  }

  if (tuning_type == "two h") {
    p0_n_combs <- length(p0_hx)*length(p0_hy)*length(p0_g)
  } else if (tuning_type == "one h") {
    p0_n_combs <- length(p0_hx)*length(p0_g)
  }

  if (is.null(p0_man_params)) {
    if (length(p0_n_combs) != 1 | !is.integer(p0_n_combs) | (p0_n_combs <= 0))
      stop("KCCA.options: argument 'p0_n_combs' should be a positive integer.")
    if (length(p0_h) == 0 | !is.numeric(p0_h) | any(p0_h <= 0))
      stop("KCCA.options: argument 'p0_h' should be a numeric vector containing double(s) greater than zero.")
    if(!is.null(p0_h_cat)) { if (length(p0_h_cat) == 0 | !is.numeric(p0_h_cat) | any(p0_h_cat <= 0))
      stop("KCCA.options: argument 'p0_h_cat' should be a numeric vector containing double(s) greater than zero.")}
    if (length(p0_hx) == 0 | !is.numeric(p0_hx) | any(p0_h <= 0))
      stop("KCCA.options: argument 'p0_h' should be a numeric vector containing double(s) greater than zero.")
    if (length(p0_hy) == 0 | !is.numeric(p0_hy) | any(p0_hy <= 0))
      stop("KCCA.options: argument 'p0_h' should be a numeric vector containing double(s) greater than zero.")
    if (length(p0_g) == 0 | !is.numeric(p0_g) | any(p0_g <= 0))
      stop("KCCA.options: argument 'p0_h' should be a numeric vector containing double(s) greater than zero.")
    if (length(p0_d) == 0 | !is.numeric(p0_d) | any(p0_d <= 0))
      stop("KCCA.options: argument 'p1_d' should be a numeric vector containing double(s) greater than zero.")
    hpar1 <- init_KCCA_one(h=p0_h,h_cat=p0_h_cat,hx=p0_hx,hy=p0_hy,g=p0_g, d=p0_d,tuning_type=tuning_type,
                             tuning_type_bandwith=tuning_type_bandwith,n_combs=n_combs)
  } else {
    hpar1 <- check_man_params_KCCA_one(p0_man_params,"0")
    p0_h <- unique(hpar1$h)
    p0_n_combs <- nrow(hpar1)/length(p0_h)
  }

  if (tuning_type_bandwith == "two h") {
  P0 <- list(
    n_combs = p0_n_combs ,
    hpar = hpar1,
    n_h = length(p0_h)*length(p0_h_cat),
    kernel_predict = p0_kernel_predict,
    kernel_predict_cat = p0_kernel_predict_cat,
    objmethod = p0_objmethod
  )
  } else if (tuning_type_bandwith == "one h") {
    P0 <- list(
      n_combs = p0_n_combs ,
      hpar = hpar1,
      n_h = length(p0_h),
      kernel_predict = p0_kernel_predict,
      kernel_predict_cat = p0_kernel_predict_cat,
      objmethod = p0_objmethod
    )
  }



  opts <- list(
    print.details = print.details,
    print.details_tuning = print.details_tuning,
    scaling = scaling,
    n_fold = n_fold,
    exactMatch = exactMatch,
    names.CV_Cat_prio = names.CV_Cat_prio,
    n_fold = n_fold,

    # Arguments for CCA or KCCA
    #d = d,
    rot = rot,
    type_predict = type_predict,

    # Argument for the tuning
    P0 = P0,
    par = par,
    par_fold = par_fold,
    nc = nc
  )
  class(opts) <- "KCCA.options.one.phase"
  return(opts)
}



init_KCCA_one <- function(h,h_cat,hx,hy,g,d,tuning_type,tuning_type_bandwith,n_combs) {

  if (tuning_type=="two h") {
    if (tuning_type_bandwith=="two h") {
      out <- merge(as.data.frame(h),as.data.frame(h_cat))
    } else if (tuning_type_bandwith=="one h") {
      out <- data.frame(h=h,h_cat=h)
    }
  out <- merge(out,as.data.frame(hx))
  out <- merge(out,as.data.frame(hy))
  out <- merge(out,as.data.frame(g))
  out <- merge(out,as.data.frame(d))
  } else if (tuning_type=="one h") {
    out <- data.frame(hx=hx,hy=hx) #Same h for x and y
    if (tuning_type_bandwith=="two h") {
      out0 <- merge(as.data.frame(h),as.data.frame(h_cat))
    } else if (tuning_type_bandwith=="one h") {
      out0 <- data.frame(h=h,h_cat=h)
    }
    out <- merge(out0,out)
    out <- merge(out,as.data.frame(g))
    out <- merge(out,as.data.frame(d))
  } else if (tuning_type=="random") {
    out <- data.frame(
      hx = sample(hx, size = n_combs, replace = TRUE),
      hy = sample(hy, size = n_combs, replace = TRUE),
      g = sample(g, size = n_combs, replace = TRUE),
      d = sample(d, size = n_combs, replace = TRUE))

    if (tuning_type_bandwith=="two h") {
      out0 <- merge(as.data.frame(h),as.data.frame(h_cat))
    } else if (tuning_type_bandwith=="one h") {
      out0 <- data.frame(h=h,h_cat=h)
    }
    out <- merge(out0, out)
  }

  names(out) <- c("h","h_cat","hx","hy","g","d")
  return(out)
}


check_man_params_KCCA_one <- function(params, phase = "0") {
  phase <- match.arg(phase)
  if (!is.data.frame(params))
    stop("KCCA.options: argument 'man_params_", phase, "' should be a data.frame object.")
  # Check if there is any missing column
  req_names <- c("h", "h_cat", "hx","hy","g", "d")
  check <- !req_names %in% colnames(params)
  missing <- paste0("'", req_names[check],"'", collapse = ", ")
  if (any(check))
    stop("KCCA.options: column(s) ", missing, " missing in 'man_params_", phase, "' data.frame.")
  # Check numeric format
  check_names <-  c("h", "h_cat", "hx","hy","g", "d")
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
      stop("number h must be the same for each quadruplet of hx, hy, g and d")
    }
  }
}
