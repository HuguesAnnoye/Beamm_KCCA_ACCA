#' @title statMatch.superOM.stabilizer
#' @description Main function to perform statistical matching using superOM
#' @param df.rec a data.frame, the recipient data set.
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @param zero.constraints a character vector with the names of the non-common variables subject to the zero constraint.
#' @param don.weights a numeric vector, the sample weights of the donor data set.
#' @param exactMatch a logical, whether or not performing exact matching (before statistical matching).
#' @param opts a list, supplied using the function \code{statMatch.superOM.options()}, containing all the options.
#' @param n_trials an integer, the number of trials used to stabilize the results.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#'
#' @importFrom dplyr intersect setdiff select mutate_all select_all all_of
#' @importFrom magrittr %>%


statMatch.superOM.stabilizer <- function(df.rec, df.don, names.CV = NULL, names.NCV = NULL,
                                         zero.constraints = NULL, don.weights = NULL,
                                         exactMatch = FALSE, opts = NULL, n_trials = 5) {
  sum_dist_comp_donor <- NULL
  pred_phase1 <- list()
  pred_phase2 <- list()
  df.phase2_avg <- 0
  phase1_tuneParams <- NULL
  phase2_tuneParams <- NULL
  for (i in seq_len(n_trials)) {
    if (opts$print.details) {
      message("\n\nSuperOM stabilizer - iteration ", i, " of ", n_trials)
    }
    tmp <- statMatch(df.rec, df.don, names.CV, names.NCV, zero.constraints, don.weights, exactMatch, opts, method = "superOM")
    sum_dist_comp_donor <- cbind(sum_dist_comp_donor, tmp$phase1$sum_dist_comp_donor)
    pred_phase1[[i]] <- tmp$phase1$df.match.don
    pred_phase2[[i]] <- tmp$phase2$df.match.don
    df.phase2_avg <- df.phase2_avg + tmp$phase2$df.match.don
    phase1_tuneParams <- rbind(phase1_tuneParams, tmp$phase1$tune[1,])
    phase2_tuneParams <- rbind(phase2_tuneParams, tmp$phase2$tune[1,])
  }

  # Phase 1
  sel_trial <- apply(sum_dist_comp_donor, 1, which.min)
  df.phase1_avg <- NULL
  for (i in seq_len(nrow(df.rec))) {
    df.phase1_avg <- rbind(df.phase1_avg, pred_phase1[[sel_trial[i]]][i,])
  }

  # Phase 2
  df.phase2_avg <- df.phase2_avg / n_trials
  ### Apply the zero constraint to obtain the final matched data set
  for (name in zero.constraints) {
    df.phase2_avg[df.phase1_avg[,paste0("ZC_", name)] == 1, name] <- 0
  }
  final_match <- df.rec %>% select(all_of(names.CV))
  final_match <- cbind(final_match, df.phase1_avg)
  final_match <- cbind(final_match, df.phase2_avg)

  res <- list(
    match = final_match,
    phase1 = list(pred = pred_phase1,
                  tune = phase1_tuneParams),
    phase2 = list(pred = pred_phase2,
                  tune = phase2_tuneParams)
  )
}
