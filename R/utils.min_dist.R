#' @title Minimal distance between two matrice.
#' @description min_dist
#' @param mtx1 a numeric matrix.
#' @param mtx2 a numeric matrix.
#' @return A numeric vect
#' @export
#'
min_dist <- function(mtx1,mtx2)
{
  res <- compute_hmin_kcca_rcpp(mtx1, mtx2, 0)

}
