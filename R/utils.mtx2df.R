#' @title mtx2df
#' @description Transform a matrix into a data.frame (dummy variables variables are transformed into factors)
#' @param mtx a numeric matrix.
#' @param var_names a character vector, for the identification of the dummy variables belonging to a given factor.
#' @return A data.frame.
#' @importFrom stringr str_detect
#' @export

mtx2df <- function(mtx, var_names) {
  df <- data.frame(matrix(nrow = nrow(mtx), ncol = 0))
  var_names_raw <- colnames(mtx)
  for (i in seq_len(length(var_names))) {
    tmp <- str_detect(var_names_raw, paste0("^", var_names[i], "$"))
    check <- sum(tmp)
    if (check == 1) {
      df <- cbind(df, mtx[, tmp])
    } else if (check == 0) {
      tmp <- str_detect(var_names_raw, paste0("^", var_names[i], "\\.\\.\\."))
      lev <- sub(".*\\.\\.\\.", "", var_names_raw[tmp])
      mtx_tmp <- kohonenBEAMM::classmat2classvec(mtx[, tmp, drop = F])
      levels(mtx_tmp) <- lev
      df <- cbind(df, mtx_tmp)
    }
  }
  colnames(df) <- var_names
  return(df)
}
