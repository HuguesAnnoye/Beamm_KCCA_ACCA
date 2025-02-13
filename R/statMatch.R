#' @title statMatch
#' @description Main function to perform statistical matching
#' @param df.rec a data.frame, the recipient data set.
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param opts a list containing all the options supplied by the function \code{statMatch.KCCA.options()}, \code{statMatch.superOM.options()}, \code{statMatch.ACCA.options()}, ..., depending on the selected method.
#' @param method a string, the method used to perform statistical matching (either \code{KCCA}, \code{superOM}, \code{ACCA}, \code{CCA}, \code{PCA} or \code{MMLP}).
#' @return The output is a list containing the following elements:
#'   \describe{
#'     \item{match}{Final matched data set (where the predicted zero constraint are applied to the non-common continuous variable(s), if any.)}
#'     \item{final}{Final matched data set (where the predicted zero constraint are applied to the non-common continuous variable(s), if any.) and where we keep the variables that are not used as common variables in the receiver}
#'     \item{phase1}{The results of phase1 of the algorthm. A list composed by:}
#'     \itemize{
#'          \item{\code{tune}, a data.frame with the combination of hyper-parameters used for cross-validation and the values of the corresponding objective function.}
#'          \item{\code{fit}, result of the fitted KCCA, CCA, ACCA, super-OM, ... with the best combination of hyper-parameters.}
#'          \item{\code{df.match.don}}
#'          \item{\code{df.match}}
#'          \item{\code{sum_dist_comp_donor}}
#'     }
#'     \item{phase2}{The results of phase2 of the algorithm. A list containing:}
#'     \itemize{
#'          \item{\code{tune}, a data.frame with the combination of hyper-parameters used for cross-validation and the values of the corresponding objective function.}
#'          \item{\code{fit}, result of the fitted KCCA, CCA, ACCA, super-OM, ... with the best combination of hyper-parameters.}
#'          \item{\code{df.match.don}}
#'          \item{\code{df.match}}
#'     }
#'     \item{opts}{options provided to the statistical matching algorithm.}
#'   }
#' @details DETAILS
#'
#' @importFrom dplyr intersect setdiff select mutate_all select_all all_of
#' @importFrom magrittr %>%
#'
#' @export

statMatch <- function(df.rec, df.don, names.CV = NULL, names.NCV = NULL,
                      zero.constraints = NULL, don.weights = NULL, opts = NULL,
                      method = c("KCCA", "superOM", "ACCA", "CCA", "PCA", "MLP")) {

  if(!is.null(opts)) opts$method <- match.arg(method)

  # Check if options exists or create default ones
  opts <- check.options(opts, method)

  if (opts$print.details) message("Statistical matching with ", method, ". Checking input arguments...")

  # Check sample weights
  don.weights <- check.don.weights(df.don, don.weights)

  # Check names & data
  check <- check.names_data(df.rec, df.don, names.CV, names.NCV, zero.constraints)
  names.CV <- check$names.CV
  names.NCV <- check$names.NCV

  # Keep the required variables only and check for NAs
  df.rec2 <- df.rec %>% select(-any_of(names.CV))
  df.rec <- df.rec %>% select(all_of(names.CV))
  df.don <- df.don %>% select(all_of(c(names.CV, names.NCV)))
  if (any(is.na(df.rec))) stop("Common variables in 'df.rec' should not contain NAs.")
  if (any(is.na(df.don))) stop("Common variables in 'df.don' should not contain NAs.")

  # Add dummy variables as zero constraints for the continuous variables
  df.don <- check.zero_constraints(df.don, zero.constraints, names.NCV)

  ### Phase 1: match the factor variables only (including zero constraints)
  if (!is.null(zero.constraints)) names.NCV_phase1 <- c(names.NCV, paste0("ZC_", zero.constraints))
  else names.NCV_phase1 <- names.NCV

  # Run: phase 1
  switch(
    opts$method,
    superOM = {
      phase1 <- statMatch.superOM.phase1(df.rec, df.don, don.weights, names.CV, names.NCV_phase1, opts)
    },
    MLP = {
      phase1 <- statMatch.MLP.phase1(df.rec, df.don, don.weights, names.CV, names.NCV_phase1, opts)
    },
    ACCA = {
      phase1 <- statMatch.ACCA.phase1(df.rec, df.don, don.weights, names.CV, names.NCV_phase1, opts)
    },
    KCCA = {
      phase1 <- statMatch.KCCA.phase1(df.rec, df.don, don.weights, names.CV, names.NCV_phase1, opts)
    },
    CCA = {
      phase1 <- statMatch.CCA.phase1(df.rec, df.don, don.weights, names.CV, names.NCV_phase1, opts)
    },
    PCA = {
      phase1 <- statMatch.PCA.phase1(df.rec, df.don, don.weights, names.CV, names.NCV_phase1, opts)
    }
  )

  ### Phase 2: match the remaining numeric variables
  names.CV_phase2 <- names(phase1$df.match)
  names.NCV_phase2 <- setdiff(names.NCV, names.CV_phase2)

  # Run: phase 2
  switch(
    opts$method,
    superOM = {
      phase2 <- statMatch.superOM.phase2(phase1$df.match, df.don, don.weights,
                                         names.CV_phase2, names.NCV_phase2, zero.constraints, opts)
    },
    MLP = {
      phase2 <- statMatch.MLP.phase2(phase1$df.match, df.don, don.weights,
                                     names.CV_phase2, names.NCV_phase2, zero.constraints, opts)
    },
    ACCA = {
      phase2 <- statMatch.ACCA.phase2(phase1$df.match, df.don, don.weights,
                                      names.CV_phase2, names.NCV_phase2, zero.constraints, opts)
    },
    KCCA = {
      phase2 <- statMatch.KCCA.phase2(phase1$df.match, df.don, don.weights,
                                      names.CV_phase2, names.NCV_phase2, zero.constraints, opts)
    },
    CCA = {
      phase2 <- statMatch.CCA.phase2(phase1$df.match, df.don, don.weights,
                                     names.CV_phase2, names.NCV_phase2, zero.constraints, opts)
    },
    PCA = {
      phase2 <- statMatch.PCA.phase2(phase1$df.match, df.don, don.weights,
                                     names.CV_phase2, names.NCV_phase2, zero.constraints, opts)
    }
  )

  ### Apply the zero constraint dummies to obtain the final matched data set
  ZC <- phase2$df.match %>% select(starts_with("ZC_"))
  match <- phase2$df.match %>% select(-starts_with("ZC_"))
  for (name in zero.constraints) {
    match[ZC[,paste0("ZC_", name)] == 1, name] <- 0
  }

  ### Output
  res <- list(match = match,
              final = cbind(df.rec2,match),
              phase1 = phase1,
              phase2 = phase2,
              opts = opts)
}
