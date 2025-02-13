#' @title scale_mtx
#' @description Scale a matrix (z-score) such that column means and standard deviations are equal to 0 and 1, respectively.
#' @param mtx a numeric matrix.
#' @param wt a numeric vector, sample weights.
#' @param center a numeric vector.
#' @param scale a numeric vector.
#' @return A numeric matrix.
#' @details If \code{center} or/and \code{scale} are equal to \code{NULL}, the function compute the column means and standard deviations.
#' @export


scale_mtx <- function(mtx, wt = NULL, center = NULL, scale = NULL) {
  n_row <- nrow(mtx)
  n_col <- ncol(mtx)
  if (is.null(center) | is.null(scale)) {
    if (!is.null(wt) & length(wt) != n_row)
      stop("length(wt) should be equal to nrow(mtx).")
    center <- wtd_col_means_rcpp(mtx, wt)
    scale <- wtd_col_sd_rcpp(mtx, wt)
    if (any(scale == 0)) scale[scale == 0] <- 1
    mtx_scal <- mtx - t(matrix(center, n_col, n_row))
    mtx_scal <- mtx_scal / t(matrix(scale, n_col, n_row))
  } else {
    if (length(center) != n_col)
      stop("length(center) should be equal to ncol(mtx).")
    if (length(scale) != n_col)
      stop("length(scale) should be equal to ncol(mtx).")
    mtx_scal <- mtx - t(matrix(center, n_col, n_row))
    mtx_scal <- mtx_scal / t(matrix(scale, n_col, n_row))
  }
  attr(mtx_scal, "scaled:center") <- center
  attr(mtx_scal, "scaled:scale") <- scale
  return(mtx_scal)
}
