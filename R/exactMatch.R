#' @title Exact Matching: SOM Tuning
#' @description What is doing this function?
#' @param df.rec a data.frame, the recipient data set.
#' @param df.don a data.frame, the donor data set.
#' @param names.CV a character vector with the names of the common variables used for statistical matching.
#' @param names.NCV a character vector with the names of the non-common variables used for statistical matching.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#'
#' @importFrom dplyr group_by_at left_join summarise_all
#' @importFrom rlang .data
#'
#' @export

exactMatch <- function(df.rec, df.don, names.CV, names.NCV){
  # Create df.don.raw and df.rec.raw, where factors are transformed into dummies
  df.rec.raw <- data.frame(df2mtx(df.rec[,names.CV]))
  df.don.raw <- data.frame(df2mtx(df.don[,c(names.CV, names.NCV)]))
  names.CV.raw <- colnames(df.rec.raw)[grepl(paste0(names.CV, "\\.\\.\\.", collapse = "|"),
                                             colnames(df.rec.raw)) |
                                         grepl(paste0(names.CV, "$", collapse = "|"),
                                               colnames(df.rec.raw))]
  # Perform exact matching and keep only the non-common variables
  message("Performing exact matching ... ")
  df.don.raw.grouped <- df.don.raw %>%
    group_by_at(names.CV.raw) %>%
    summarise_all(~mean(.data, na.rm = TRUE))
  df.don.matched.raw <- left_join(x = df.rec.raw,
                                  y = df.don.raw.grouped,
                                  by = names.CV.raw) %>%
    select(-which(names(.data) %in% names.CV.raw))
  n_exact_match <- sum(rowSums(is.na(df.don.matched.raw)) < ncol(df.don.matched.raw))
  message("   Done ! In the recevier dataset, ",n_exact_match, " observations out of ",nrow(df.rec), " had an exact match with the donor dataset.")
  df.don.matched <- mtx2df(as.matrix(df.don.matched.raw), var_names = names.NCV)
  return(df.don.matched)
}
