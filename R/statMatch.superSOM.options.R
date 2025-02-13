#' @title statMatch.superOM.options
#' @description Create a list with the options for the statistical matching algorithm based on superOM
#' @param print.details a logical, whether or not details about the different procedures involved in the statistical matching are printed on the terminal.
#' @param scaling a character string, the name of the method used to scale the variables in the receiver and donor data sets: \code{"z-score"}, \code{"min-max"} or \code{"no"} (no scaling).
#' @param tuning_type a character string, the method used to find the optimal tuning parameters: either \code{"grid"} or \code{"random"} search.
#' @param p1_rand_combs a positive integer, number of random combinations of tuning parameters generated in phase 1 (used only when \code{tuning_type = "random"}).
#' @param p2_rand_combs a positive integer, number of random combinations of tuning parameters generated in phase 2 (used only when \code{tuning_type = "random"}).
#' @param p1_square_grid a logical, whether or not a superOM with a square grid is used in phase 1.
#' @param p2_square_grid a logical, whether or not a superOM with a square grid is used in phase 2.
#' @param p1_dim_min a positive integer, minimal length of the x and y axis of the superOM grid in phase 1.
#' @param p2_dim_min a positive integer, minimal length of the x and y axis of the superOM grid in phase 2.
#' @param p1_dim_max a positive integer, maximal length of the x and y axis of the superOM grid in phase 1.
#' @param p2_dim_max a positive integer, maximal length of the x and y axis of the superOM grid in phase 2.
#' @param p1_dim_step a positive integer, step size of the x and y axis of the superOM grid in phase 1 (used only when \code{tuning_type = "grid"}).
#' @param p2_dim_step a positive integer, step size of the x and y axis of the superOM grid in phase 2 (used only when \code{tuning_type = "grid"}).
#' @param p1_wCV_min a double, minimal value of the weight of the superOM layer containing the common variables in phase 1.
#' @param p2_wCV_min a double, minimal value of the weight of the superOM layer containing the common variables in phase 2.
#' @param p1_wCV_max a double, maximal value of the weight of the superOM layer containing the common variables in phase 1.
#' @param p2_wCV_max a double, maximal value of the weight of the superOM layer containing the common variables in phase 2.
#' @param p1_wCV_step a double, step size of the weight of the superOM layer containing the common variables in phase 1 (used only when \code{tuning_type = "grid"}).
#' @param p2_wCV_step a double, step size of the weight of the superOM layer containing the common variables in phase 2 (used only when \code{tuning_type = "grid"}).
#' @param n_fold a positive integer, number of folds used for the cross-validation of the tuning parameters.
#' @param n_nodes a positive integer, number of parallel sessions (CPUs) used for the cross-validation of the tuning parameters.
#' @param topo a character string, topology of the superOM grid, either \code{"hexagonal"} or \code{"rectangular"}.
#' @param neighbourhood.fct a character string, neighbourhood function used to train the superOM, either \code{"gaussian"} or \code{"bubble"}.
#' @param toroidal a logical, whether the superOM grid is toroidal or not.
#' @param rlen a positive integer, the number of times the complete data set will be presented to the network when the training of the superOM.
#' @param dist.fcts a character string, distance function used for the individual data layers (only \code{"euclidean"}).
#' @param maxNA.fraction a double, the maximal fraction of values that may be NA to prevent the row to be removed when training the superOM.
#' @param seed_phase1 a positive integer, the seed for random number generation in phase 1 (a random value if \code{NULL}).
#' @param seed_phase2 a positive integer, the seed for random number generation in phase 2 (a random value if \code{NULL}).
#' @param man_params_1 a data.frame, a manual way to specify the grid of tuning parameters in phase 1.
#' @param man_params_2 a data.frame, a manual way to specify the grid of tuning parameters in phase 2.
#' @return A list with the options required to the statistical matching algorithm using superOM.
#'
#' @export

statMatch.superOM.options <- function(print.details = TRUE,
                                      scaling = c("z-score", "min-max", "no"),
                                      # Arguments for tuning
                                      tuning_type = c("grid", "random"),
                                      p1_rand_combs = 10L,
                                      p2_rand_combs = 10L,
                                      p1_square_grid = TRUE,
                                      p2_square_grid = TRUE,
                                      p1_dim_min = 2L,
                                      p2_dim_min = 2L,
                                      p1_dim_max = 20L,
                                      p2_dim_max = 20L,
                                      p1_dim_step = 1L,
                                      p2_dim_step = 1L,
                                      p1_wCV_min = 0,
                                      p2_wCV_min = 0,
                                      p1_wCV_max = 1,
                                      p2_wCV_max = 1,
                                      p1_wCV_step = 0.1,
                                      p2_wCV_step = 0.1,
                                      n_fold = 5L,
                                      n_nodes = 1L,
                                      # Arguments used by kohonen
                                      topo = c("hexagonal", "rectangular"),
                                      neighbourhood.fct = c("gaussian", "bubble"),
                                      toroidal = FALSE,
                                      rlen = 250L,
                                      dist.fcts = c("euclidean"),
                                      maxNA.fraction = 0.9,
                                      seed_phase1 = NULL,
                                      seed_phase2 = NULL,
                                      man_params_1 = NULL,
                                      man_params_2 = NULL
                                      ) {

  scaling <- match.arg(scaling)
  tuning_type <- match.arg(tuning_type)
  topo <- match.arg(topo)
  neighbourhood.fct <- match.arg(neighbourhood.fct)
  dist.fcts = match.arg(dist.fcts)

  # Check inputs

  if (length(print.details) != 1 | !is.logical(print.details) | any(is.na(print.details)))
    stop("superOM.options: argument 'print.details' should be a logical.")
  if (length(n_fold) != 1 | !is.integer(n_fold) | (n_fold <= 0))
    stop("superOM.options: argument 'n_fold' should be a positive integer.")
  if (length(n_nodes) != 1 | !is.integer(n_nodes) | (n_nodes <= 0))
    stop("superOM.options: argument 'n_nodes' should be a positive integer.")
  if (length(toroidal) != 1 | !is.logical(toroidal) | any(is.na(toroidal)))
    stop("superOM.options: argument 'toroidal' should be a logical.")
  if (length(rlen) != 1 | !is.integer(rlen) | (rlen <= 0))
    stop("superOM.options: argument 'rlen' should be a positive integer.")
  if (length(maxNA.fraction) != 1 | !is.double(maxNA.fraction) | (maxNA.fraction < 0) | (maxNA.fraction > 1))
    stop("superOM.options: argument 'maxNA.fraction' should be a double between 0 and 1.")

  if (is.null(seed_phase1)) seed_phase1 <- sample(1:1e9, 1)
  else if (length(seed_phase1) != 1 | !is.integer(seed_phase1) | (seed_phase1 <= 0))
    stop("superOM.options: argument 'seed_phase1' should be a positive integer.")
  if (is.null(seed_phase2)) seed_phase2 <- sample(1:1e9, 1)
  else if (length(seed_phase2) != 1 | !is.integer(seed_phase2) | (seed_phase2 <= 0))
    stop("superOM.options: argument 'seed_phase2' should be a positive integer.")

  if (is.null(man_params_1)) {
    if (length(p1_square_grid) != 1 | !is.logical(p1_square_grid) | any(is.na(p1_square_grid)))
      stop("superOM.options: argument 'p1_square_grid' should be a logical.")
    if (length(p1_rand_combs) != 1 | !is.integer(p1_rand_combs) | (p1_rand_combs <= 0))
      stop("superOM.options: argument 'p1_rand_combs' should be a positive integer.")
    if (length(p1_dim_min) != 1 | !is.integer(p1_dim_min) | (p1_dim_min <= 0))
      stop("superOM.options: argument 'p1_dim_min' should be a positive integer.")
    if (length(p1_dim_max) != 1 | !is.integer(p1_dim_max) | (p1_dim_max <= 0) | (p1_dim_max < p1_dim_min))
      stop("superOM.options: argument 'p1_dim_max' should be a positive integer, higher than or equal to 'p1_xdim_min'.")
    if (length(p1_dim_step) != 1 | !is.integer(p1_dim_step) | (p1_dim_step <= 0))
      stop("superOM.options: argument 'p1_dim_step' should be a positive integer.")
    if (length(p1_wCV_min) != 1 | !is.numeric(p1_wCV_min) | (p1_wCV_min < 0) | (p1_wCV_min > 1))
      stop("superOM.options: argument 'p1_wCV_min' should be a double between 0 and 1.")
    if (length(p1_wCV_max) != 1 | !is.numeric(p1_wCV_max) | (p1_wCV_max < 0) | (p1_wCV_max > 1) | (p1_wCV_max < p1_wCV_min))
      stop("superOM.options: argument 'p1_wCV_max' should be a double between 0 and 1, higher than or equal to 'p1_wCV_min'.")
    if (length(p1_wCV_step) != 1 | !is.numeric(p1_wCV_step) | (p1_wCV_step < 0) | (p1_wCV_step > 1))
      stop("superOM.options: argument 'p1_wCV_step' should be a double between 0 and 1.")
    if (tuning_type == "grid") {
      p1_params <- init_grid_search_superOM(p1_dim_min, p1_dim_max, p1_dim_step, p1_square_grid,
                                            p1_wCV_min, p1_wCV_max, p1_wCV_step)
    } else {
      p1_params <- init_random_search_superOM(p1_dim_min, p1_dim_max, p1_square_grid,
                                              p1_wCV_min, p1_wCV_max, p1_rand_combs)
    }
  } else {
    p1_params <- check_man_params_superOM(man_params_1, phase = "1")
  }

  if (is.null(man_params_2)) {
    if (length(p2_square_grid) != 1 | !is.logical(p2_square_grid) | any(is.na(p2_square_grid)))
      stop("superOM.options: argument 'p2_square_grid' should be a logical.")
    if (length(p2_rand_combs) != 1 | !is.integer(p2_rand_combs) | (p2_rand_combs <= 0))
      stop("superOM.options: argument 'p2_rand_combs' should be a positive integer.")
    if (length(p2_dim_min) != 1 | !is.integer(p2_dim_min) | (p2_dim_min <= 0))
      stop("superOM.options: argument 'p2_dim_min' should be a positive integer.")
    if (length(p2_dim_max) != 1 | !is.integer(p2_dim_max) | (p2_dim_max <= 0) | (p2_dim_max < p2_dim_min))
      stop("superOM.options: argument 'p2_dim_max' should be a positive integer, higher than or equal to 'p2_xdim_min'.")
    if (length(p2_dim_step) != 1 | !is.integer(p2_dim_step) | (p2_dim_step <= 0))
      stop("superOM.options: argument 'p2_dim_step' should be a positive integer.")
    if (length(p2_wCV_min) != 1 | !is.numeric(p2_wCV_min) | (p2_wCV_min < 0) | (p2_wCV_min > 1))
      stop("superOM.options: argument 'p2_wCV_min' should be a double between 0 and 1.")
    if (length(p2_wCV_max) != 1 | !is.numeric(p2_wCV_max) | (p2_wCV_max < 0) | (p2_wCV_max > 1) | (p2_wCV_max < p2_wCV_min))
      stop("superOM.options: argument 'p2_wCV_max' should be a double between 0 and 1, higher than or equal to 'p2_wCV_min'.")
    if (length(p2_wCV_step) != 1 | !is.numeric(p2_wCV_step) | (p2_wCV_step < 0) | (p2_wCV_step > 1))
      stop("superOM.options: argument 'p2_wCV_step' should be a double between 0 and 1.")
    if (tuning_type == "grid") {
      p2_params <- init_grid_search_superOM(p2_dim_min, p2_dim_max, p2_dim_step, p2_square_grid,
                                            p2_wCV_min, p2_wCV_max, p2_wCV_step)
    } else {
      p2_params <- init_random_search_superOM(p2_dim_min, p2_dim_max, p2_square_grid,
                                              p2_wCV_min, p2_wCV_max, p2_rand_combs)
    }
  } else {
    p2_params <- check_man_params_superOM(man_params_2, phase = "2")
  }

  opts <- list(
    print.details = print.details,
    scaling = scaling,
    tuning_type = tuning_type,
    p1_params = p1_params,
    p2_params = p2_params,
    n_fold = n_fold,
    n_nodes = n_nodes,
    topo = topo,
    neighbourhood.fct = neighbourhood.fct,
    toroidal = toroidal,
    rlen = rlen,
    dist.fcts = dist.fcts,
    maxNA.fraction = maxNA.fraction,
    seed_phase1 = seed_phase1,
    seed_phase2 = seed_phase2
  )

  class(opts) <- "superOM.options"

  return(opts)
}


init_grid_search_superOM <- function(dim_min, dim_max, dim_step, is_square, wCV_min, wCV_max, wCV_step) {
  if (is_square) {
    out <- expand.grid(
      xdim = seq(dim_min, dim_max, by = dim_step),
      weight.CV = seq(wCV_min, wCV_max, by = wCV_step)
    )
    out$ydim <- out$xdim
  } else {
    out <- expand.grid(
      xdim = seq(dim_min, dim_max, by = dim_step),
      ydim = seq(dim_min, dim_max, by = dim_step),
      weight.CV = seq(wCV_min, wCV_max, by = wCV_step)
    )
  }
  return(out)
}

init_random_search_superOM <- function(dim_min, dim_max, is_square, wCV_min, wCV_max, n_combs) {
  if (is_square) {
    out <- data.frame(
      xdim = sample(dim_min:dim_max, size = n_combs, replace = T),
      weight.CV = runif(n = n_combs, min = wCV_min, max = wCV_max)
    )
    out$ydim <- out$xdim
  } else {
    out <- data.frame(
      xdim = sample(dim_min:dim_max, size = n_combs, replace = T),
      ydim = sample(dim_min:dim_max, size = n_combs, replace = T),
      weight.CV = runif(n = n_combs, min = wCV_min, max = wCV_max)
    )
  }
  return(out)
}

check_man_params_superOM <- function(params, phase = c("1", "2")) {
  phase <- match.arg(phase)
  if (!is.data.frame(params)) stop("superOM.options: argument 'man_params_", phase, "' should be a data.frame object.")
  req_names <- c("xdim", "ydim", "weight.CV")
  check <- !req_names %in% colnames(params)
  missing <- paste0("'", req_names[check],"'", collapse = ", ")
  if (any(check))
    stop("superOM.options: column(s) ", missing, " missing in 'man_params_", phase, "' data.frame.")
  if (!is.integer(params$xdim) | (params$xdim <= 0L))
    stop("superOM.options: column 'xdim' in 'man_params_", phase, "' data.frame should contain positive integer(s).")
  if (!is.integer(params$ydim) | (params$ydim <= 0L))
    stop("superOM.options: column 'ydim' in 'man_params_", phase, "' data.frame should contain positive integer(s).")
  if (!is.numeric(params$weight.CV) | (params$weight.CV <= 0) | (params$weight.CV >= 1))
    stop("superOM.options: column 'weight.CV' in 'man_params_", phase, "' data.frame should contain double(s) between 0 and 1 (not included).")
  return(params)
}
