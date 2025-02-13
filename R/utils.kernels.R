#' @title Contain all the kernels functions
#' @details DETAILS
#' @param x first (or unique) element in the kernel
#' @param y second element in the kernel
#' @param h a numerical value, the bandwidth
#' @importFrom stats contrasts model.matrix
#' @noRd


gausskcca <- function(x, y, h) {
  (1 / sqrt(2 * pi)) * exp(-((x - y) %*% (x - y)) / (2 * h^2))
}
Kgauss <- function(x, y, h) {
  t(apply(x, 1, function(p) apply(y, 1, function(q) gausskcca(p, q, h))))
}

gausskcca_new <- function(x) {
  (2 * pi)^(-1 / 2) * (exp((-x^2) / (2)))
}
epan_new <- function(x) {
  (3 / 4) * (1 - (x)^2) * ind(((x)^2), 1, 0)
}

gauss2d <- function(x, h) (1 / sqrt(2 * pi)) * exp(-((x[[1]] - x[[2]]) / h)^2)
ind <- function(x, a, b) {
  ifelse(b <= x & x <= a, 1, 0)
}
epan2d <- function(x, h) {
  (3 / 4) * (1 - ((x[[1]] - x[[2]]) / h)^2) * ind((((x[[1]] - x[[2]]) / h)^2), 1, 0)
}
unif2d <- function(x, h) {
  1 / 2 * ind((((x[[1]] - x[[2]]) / h)^2), 1, 0)
}

gauss3d <- function(x) ((2 * pi)^(-1) * (exp((-x[[1]]^2) / (2))) * (exp((-x[[2]]^2) / (2))))

epan3d <- function(x) ((9 / 16) * (1 - x[[1]]^2) * ((x[[1]]^2) < 1) * (1 - x[[2]]^2) * ((x[[2]]^2) < 1))
# epan3d <- function(x) ((9 / 16) * (1 - x[[1]]^2) * ((x[[1]]^2) < 1) * (1 - x[[2]]^2) * ((x[[2]]^2) < 1))

unif3d <- function(x) {
  1 / 4 * ind(((x[[1]])^2), 1, 0) * ind(((x[[2]])^2), 1, 0)
}

prepare_data_cca <- function(df1, df2, weights = NULL, dfs = NULL) {
  if (is.null(weights)) {
    nomx <- colnames(df1)
    nomy <- colnames(df2)
    df <- cbind(df1, df2)
    df <- df[apply(df, 1, function(x) !any(is.na(x))), ]
    df3 <- df[, nomx]
    df4 <- df[, nomy]
    # current.na.action <- options('na.action')
    # options(na.action='na.pass')
    df5 <- model.matrix(~., df3)
    df6 <- model.matrix(~., df4)
    df5 <- df5[, -1]
    df6 <- df6[, -1]

    df7 <- model.matrix(~., df3, contrasts.arg = lapply(data.frame(df3[, sapply(data.frame(df3), is.factor)]), contrasts, contrasts = FALSE))
    df8 <- model.matrix(~., df4, contrasts.arg = lapply(data.frame(df4[, sapply(data.frame(df4), is.factor)]), contrasts, contrasts = FALSE))
    df7 <- df7[, -1]
    df8 <- df8[, -1]
    # options('na.action' = current.na.action$na.action)
    if (is.null(dfs)) {
      return(list(mtx.don.CV.raw = df5, mtx.don.NCV.raw = df6, mtx.don.weights.raw = NULL, mtx.don.CV.all = df7, mtx.don.NCV.all = df8))
    } else {
      df9 <- model.matrix(~., dfs)
      df9 <- df9[, -1]
      return(list(mtx.don.CV.raw = df5, mtx.don.NCV.raw = df6, mtx.don.weights.raw = NULL, mtx.don.CV.all = df7, mtx.don.NCV.all = df8, mtx.rec.CV.raw = df9))
    }
  }
  else {
    nomx <- colnames(df1)
    nomy <- colnames(df2)
    nomw <- colnames(weights)
    df <- cbind(df1, df2, weights)
    df <- df[apply(df, 1, function(x) !any(is.na(x))), ]
    df3 <- df[, nomx]
    df4 <- df[, nomy]
    weightsmm <- df[, nomw]
    df5 <- model.matrix(~., df3)
    df6 <- model.matrix(~., df4)
    df5 <- df5[, -1]
    df6 <- df6[, -1]
    df7 <- model.matrix(~., df3, contrasts.arg = lapply(data.frame(df3[, sapply(data.frame(df3), is.factor)]), contrasts, contrasts = FALSE))
    df8 <- model.matrix(~., df4, contrasts.arg = lapply(data.frame(df4[, sapply(data.frame(df4), is.factor)]), contrasts, contrasts = FALSE))
    df7 <- df7[, -1]
    df8 <- df8[, -1]

    if (is.null(dfs)) {
      return(list(mtx.don.CV.raw = df5, mtx.don.NCV.raw = df6, mtx.don.weights.raw = weightsmm, mtx.don.CV.all = df7, mtx.don.NCV.all = df8))
    } else {
      df9 <- model.matrix(~., dfs)
      df9 <- df9[, -1]
      return(list(mtx.don.CV.raw = df5, mtx.don.NCV.raw = df6, mtx.don.weights.raw = weightsmm, mtx.don.CV.all = df7, mtx.don.NCV.all = df8, mtx.rec.CV.raw = df9))
    }
  }
}
