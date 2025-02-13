#' @title df2mtx
#' @description Transform a data.frame into a matrix (factors are transformed into \eqn{n} dummy variables, where \eqn{n} is the number of levels)
#' @param df a data.frame.
#' @param n_1_levels a logical value, if TRUE keep only \eqn{n-1} dummies for each factor.
#' @return A numeric matrix.
#'
#' @export

df2mtx <- function(df, n_1_levels = FALSE) {
  out <- NULL
  var_names <- colnames(df)
  for (i in seq_len(ncol(df))) {
    if (is.factor(df[, i, drop = T])) {
      tmp <- kohonenBEAMM::classvec2classmat(df[, i, drop = T])
      colnames(tmp) <- paste(var_names[i], levels(df[, i, drop = T]), sep = "...")
      if (n_1_levels) {
        out <- cbind(out, tmp[, -1, drop = F])
      } else {
        out <- cbind(out, tmp)
      }
    } else {
      tmp <- df[, i, drop = T]
      out <- cbind(out, tmp)
      colnames(out)[ncol(out)] <- var_names[i]
    }
  }
  return(out)
}
