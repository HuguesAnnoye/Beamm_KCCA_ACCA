#' @title normalize_mtx_inv
#' @description Transform a matrix from scaled (min-max) to unscaled.
#' @param mtx_norm a matrix object.
#' @param min a numeric vector.
#' @param max a numeric vector.
#' @return A matrix object.
#' @details DETAILS
#' @export

normalize_mtx_inv <- function(mtx_norm, min, max) {
  n_row <- nrow(mtx_norm)
  n_col <- ncol(mtx_norm)
  mtx <- mtx_norm * t(matrix(max - min, n_col, n_row))
  mtx <- mtx + t(matrix(min, n_col, n_row))
  return(mtx)
}
