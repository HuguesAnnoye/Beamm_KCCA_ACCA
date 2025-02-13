#' @title normalize_mtx
#' @description Scale a matrix (min-max) such that column minimums and maximums are equal to 0 and 1, respectively.
#' @param mtx a numeric matrix.
#' @param min a numeric vector.
#' @param max a numeric vector.
#' @return A numeric matrix.
#' @export


normalize_mtx <- function(mtx, min = NULL, max = NULL) {
  n_row <- nrow(mtx)
  n_col <- ncol(mtx)
  if (is.null(min) | is.null(max)) {
    min <- apply(mtx, 2, min)
    max <- apply(mtx, 2, max)
    mtx_norm <- mtx - t(matrix(min, n_col, n_row))
    mtx_norm <- mtx_norm / t(matrix(max - min, n_col, n_row))
  } else {
    mtx_norm <- mtx - t(matrix(min, n_col, n_row))
    mtx_norm <- mtx_norm / t(matrix(max - min, n_col, n_row))
  }
  attr(mtx_norm, "min") <- min
  attr(mtx_norm, "max") <- max
  return(mtx_norm)
}
