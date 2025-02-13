compute_MCR <- function(y, y_hat, weights = NULL) {
  if (is.null(weights)) {
    MCRw <- 1 - mean(y == y_hat, na.rm = T)
  }
  else {
    MCRw <- 1 - wtd_mean_rcpp(y == y_hat, wt = weights)
  }
  return(MCRw)
}
