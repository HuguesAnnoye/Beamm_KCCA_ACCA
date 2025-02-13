#' @title statMatch.superOM.phase1.tuning
#' @description Tuning (phase1) of the statistical matching algorithm using superOM
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param params a data.frame, each row correspond to a different combination of hyper-parameters.
#' @param opts a list, supplied using the function \code{statMatch.superOM.options()}, containing all the options.
#' @return A data.frame, ...
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom dplyr arrange
#' @importFrom rlang .data
#' @noRd

statMatch.superOM.phase1.tuning <- function(df.don, names.CV, names.NCV, don.weights, params, opts) {

  names.CV.cat <- names.CV[unlist(lapply(df.don %>% select(all_of(names.CV)), is.factor))]

  # Prepare data for cross-validation
  input <- prepare_data_tuning(df.don, names.CV, names.NCV, don.weights, opts,
                               comp.mtx = T, names.CV.cat = names.CV.cat)

  # Initialize cluster
  nnodes <- pmin(opts$n_nodes, parallel::detectCores() - 1, nrow(params))
  cl <- parallel::makeCluster(nnodes, outfile = "")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  message("\nPhase 1. Tuning grid size and layer weights...")
  i <- NULL
  x <- foreach::foreach(i = seq_len(nrow(params)), .combine = rbind) %dopar% {
    t_start <- Sys.time()
    set.seed(opts$seed_phase1)
    PRED.raw <- NULL
    for (fold in 1:opts$n_fold) {
      tryCatch({
        fit <- statMatch.superOM.fit(input$raw_data[[fold]], params$xdim[i],
                                     params$ydim[i], params$weight.CV[i],
                                     input$weights[[fold]], opts)
        },
        error = function(e) {
          print(paste0("Error in phase 1! Problem in the kohonen package. ",
                       "grid_xdim = ", params$xdim[i],
                       ", grid_ydim = ", params$ydim[i],
                       ", weight.CV = ", params$weight.CV[i]))
          print(e)
        })
      pred <- statMatch.superOM.phase1.predict(fit,
                                               input$raw_data[[fold]],
                                               input$weights[[fold]],
                                               input$comp.idx[[fold]],
                                               opts$scaling)
      PRED.raw <- rbind(PRED.raw, pred)
    }
    val.objfun <- compute_obj_function(input$mtx.don.NCV.raw,
                                       PRED.raw,
                                       don.weights,
                                       type = "wRMSE")
    timing <- as.numeric(difftime(Sys.time(), t_start, units = 'mins'))
    # Print details of the current iteration
    if (opts$print.details)
      message("grid_xdim = ", params$xdim[i],
              ", grid_ydim = ", params$ydim[i],
              ", weight.CV = ", params$weight.CV[i],
              ", val_objfun = ", round(val.objfun, 3),
              ", n_NAs = ", sum(rowSums(is.na(PRED.raw)) > 0),
              ", time = ", timing)
    res <- c(params$xdim[i], params$ydim[i], params$weight.CV[i],
             val.objfun, sum(rowSums(is.na(PRED.raw)) > 0), timing)
  }

  if (is.null(dim(x))) {
    x <- as.data.frame(t(x))
  } else {
    x <- as.data.frame(x)
  }
  names(x) <- c("xdim", "ydim", "weight.CV", "val_objfun", "n_NAs", "time")
  rownames(x) <- NULL
  x <- x %>% arrange(.data$val_objfun)

  return(x)
}
