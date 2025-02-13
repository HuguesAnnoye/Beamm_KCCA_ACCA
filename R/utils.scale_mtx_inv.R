#' @title scale_mtx_inv
#' @description Transform a matrix from scaled (z-score) to unscaled.
#' @param mtx_scaled a numeric matrix.
#' @param center a numeric vector.
#' @param scale a numeric vector.
#' @return A numeric matrix.
#' @export

scale_mtx_inv <- function(mtx_scaled, center, scale) {
  n_row <- nrow(mtx_scaled)
  n_col <- ncol(mtx_scaled)
  mtx <- mtx_scaled * t(matrix(scale, n_col, n_row))
  mtx <- mtx + t(matrix(center, n_col, n_row))
  return(mtx)
}
