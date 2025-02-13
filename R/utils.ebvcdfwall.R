#' @title cramer von minsen discrete
#' @description Compute the value of the compute_pseudo_Rsquare ans make the mean or the sum
#' @param x Dataframe
#' @param ... ellipsis
#' @param wt a numeric vector, sample weights.
#' @param data data
#' @return a function
#' @importFrom bivariate ebvcdf
#' @export

ebvcdfwall <- function(x, wt, ..., data) {
  if (missing(data)) {
    variable.names <- NULL
    data <- NULL
    for (i in 1:ncol(x)){
    # {
      yname <- as.character(substitute(x[,i]))[1]
      y <- as.numeric(x[,i])
      variable.names  <- c(variable.names,yname)
      data <- cbind(data,y)
    }
    #data <- cbind(x, y, z)
    colnames(data) <- variable.names
  }
  else {
    data <- as.matrix(data)
    variable.names <- colnames(data)
  }

  sf <- function(x, wt) {
    sf <- .THIS()
    .cbv.evalall(sf, .ebvcdf.eval.ext.wtall, x)
  }

  new("EBVCDFWALL", sf,
      .class.info = .get.infow("EBVCDFWALL"),
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
setClass ("EBVCDFWALL", contains = c ("EBVW", "SBV", "BV", "CDF", "function") )

.ebvcdf.eval.ext.wtall <- function(sf, x) {
  n0 <- sf@n
  x0 <- sf@data
  w0 <- sf@wt

  n <- nrow(x)
  z1 <- numeric(n)
  #print(dim(x0))
  #print(dim(x))
  for (i in seq_len(n))
  {
    Ix <- (x0[,1] <= x[i,1])
    if (ncol(x0)>=2){
    for (j in 2:ncol(x0)) {
      Ix<- Ix*(x0[,j] <= x[i,j])
      #Ix<- Ix*(x0[,i] <= x[i,i])
   # Ix <- (x0 <= x [i])
    #Iy <- (y0 <= y [i])
    #Iz <- (z0 <= z [i])
    #
    }
    }
    z1 [i] <- sum(w0* Ix) / n0
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
  "EBVCDFW3",  "Trivariate Empirical Distribution Function with weight",         "Fh (x, y, z)",
  "EBVCDFWALL",  "Multivariate Empirical Distribution Function with weight",         "Fh (x)"
  ),, 3, byrow=TRUE)

.cbv.evalall = function (sf, evalf, x)
{
#x = as.numeric (x)
#y = as.numeric (y)
#z = as.numeric (z)
  v <- NULL
  for (i in 1:ncol(x)){
    # {
    y <- as.numeric(x[,i])
    v = cbind (v,y, deparse.level=0)
  }
#v = cbind (x, y, z, deparse.level=0)
.val.finiteall (v)
evalf (sf, v)
}

.val.finiteall = function (x, n=ncol (x))
{	if (nrow (x) == 0 || n != ncol (x) )
  stop ("\nunsuitable evaluation bins/points\n(possibly NULL or zero-length vectors)")
  if (! all (is.finite (x) ) )
    stop ("evaluation bins/points need to be finite")
}
