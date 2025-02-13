#' @title cramer von minsen discrete
#' @description Compute the value of the compute_pseudo_Rsquare ans make the mean or the sum
#' @param x Variable x
#' @param y Variable y
#' @param z Variable z
#' @param ... ellipsis
#' @param wt a numeric vector, sample weights.
#' @param data data
#' @return a function
#' @importFrom bivariate ebvcdf
#' @export

ebvcdfw3 <- function(x, y, z, wt, ..., data) {
  if (missing(data)) {
    xname <- as.character(substitute(x))[1]
    yname <- as.character(substitute(y))[1]
    zname <- as.character(substitute(z))[1]
    variable.names <- c(xname, yname, zname)
    x <- as.numeric(x)
    y <- as.numeric(y)
    y <- as.numeric(z)
    if (length(x) != length(y)) {
      stop("lengths of x and y, not same")
    }
    if (length(x) != length(z)) {
      stop("lengths of x and z, not same")
    }
    data <- cbind(x, y, z)
    colnames(data) <- variable.names
  }
  else {
    data <- as.matrix(data)
    variable.names <- colnames(data)
  }

  sf <- function(x, y, z, wt) {
    sf <- .THIS()
    .cbv.eval3(sf, .ebvcdf.eval.ext.wt3, x, y, z)
  }

  new("EBVCDFW3", sf,
      .class.info = .get.infow("EBVCDFW3"),
      variable.names = variable.names,
      #n = nrow(data),
      n  = sum(wt),
      wt = wt,
      data = data
  )
}

setClass ("EBVW", contains="SBV",
          slots = list (
            variable.names="character",
            n="numeric",
            wt="vector",
            data="matrix") )
setClass ("EBVCDFW3", contains = c ("EBVW", "SBV", "BV", "CDF", "function") )

.ebvcdf.eval.ext.wt3 <- function(sf, x, y, z) {
  n0 <- sf@n
  x0 <- sf@data [, 1]
  y0 <- sf@data [, 2]
  z0 <- sf@data [, 3]
  w0 <- sf@wt

  n <- length(x)
  z1 <- numeric(n)
  for (i in seq_len(n))
  {
    Ix <- (x0 <= x [i])
    Iy <- (y0 <= y [i])
    Iz <- (z0 <= z [i])
    z1 [i] <- sum(w0* Ix * Iy * Iz) / n0
  }
  z1
}

.get.infow = function (cn)
{	I = match (cn, .class.infow [,1])
if (is.na (I [1]) )
  stop ("constructor error")
.class.infow [I, 2:3]
}

.THIS = function ()
  sys.function (-1)

.class.infow = matrix (c (
  "DUBVPMF", "Bivariate Discrete Uniform Mass Function",           "f (x, y)",
  "BNBVPMF", "Bivariate Binomial Mass Function",                   "f (x, y)",
  "PBVPMF",  "Bivariate Poisson Mass Function",                    "f (x, y)",
  "GBVPMF",  "Bivariate Categorical Mass Function",                "f (x, y)",
  "CUBVPDF", "Bivariate Continuous Uniform Density Function",      "f (x, y)",
  "NBVPDF",  "Bivariate Normal Density Function",                  "f (x, y)",
  "BMBVPDF", "Bivariate Bimodal Density Function",                 "f (x, y)",
  "DUBVCDF", "Bivariate Discrete Uniform Distribution Function",   "F (x, y)",
  "BNBVCDF", "Bivariate Binomial Distribution Function",           "F (x, y)",
  "PBVCDF",  "Bivariate Poisson Distribution Function",            "F (x, y)",
  "CUBVCDF", "Bivariate Continuous Uniform Distribution Function", "F (x, y)",
  "NBVCDF",  "Bivariate Normal Distribution Function",             "F (x, y)",
  "BMBVCDF", "Bivariate Bimodal Distribution Function",            "F (x, y)",
  "NTVPDF",  "Trivariate Normal Density Function",                 "f (x, y, z)",
  "DTVPDF",  "Trivariate Dirichlet Density Function",              "f (x, y, z = 1 - x - y)",
  "KBVPDF",  "Bivariate Kernel Density Estimate",   "equivalent to fh (x, y)",
  "EBVCDF",  "Bivariate Empirical Distribution Function",         "Fh (x, y)",
  "EBVCDFW",  "Bivariate Empirical Distribution Function with weight",         "Fh (x, y)",
  "EBVCDFW3",  "Trivariate Empirical Distribution Function with weight",         "Fh (x, y, z)"
  ),, 3, byrow=TRUE)

.cbv.eval3 = function (sf, evalf, x, y, z)
{	x = as.numeric (x)
y = as.numeric (y)
z = as.numeric (z)
v = cbind (x, y, z, deparse.level=0)
.val.finite3 (v)
evalf (sf, v [,1], v [,2], v [,3])
}

.val.finite3 = function (x, n=3)
{	if (nrow (x) == 0 || n != ncol (x) )
  stop ("\nunsuitable evaluation bins/points\n(possibly NULL or zero-length vectors)")
  if (! all (is.finite (x) ) )
    stop ("evaluation bins/points need to be finite")
}
