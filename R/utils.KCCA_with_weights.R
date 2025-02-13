#' @title  KCCA with weights
#' @param x matrix 1.
#' @param y matrix 2.
#' @param w a vector of weights
#' @param kernel a charachter string, with the type of kernel
#' @param kpar list with the bandwidth for x
#' @param kpar2 list with the bandwidth for y
#' @param gamma regularization parameter
#' @param ... other options
#' @param ncomps number of components
#' @rawNamespace import(kernlab, except = predict)
#' @importFrom methods is new
#' @importFrom RSpectra eigs
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @noRd

setGeneric("kccaw",function(x, y, w = NULL, kernel = "rbfdot",  kpar = list(sigma = 0.1),kpar2 = list(sigma = 0.1), gamma = 0.1, ncomps = 10, ...) standardGeneric("kccaw"))

#' @rdname kccaw
#' @noRd
setMethod("kccaw", signature(x = "matrix"),
          function(x, y, w = NULL, kernel = "rbfdot", kpar = list(sigma = 0.1), kpar2 = list(sigma = 0.1), gamma = 0.1, ncomps = 10, ...)
          {
            ## mfunc is a helper to compute matrix functions
            mfunc <- function(a, fn = sqrt) {
              e <- eigen(a); y <- e$vectors; v <- e$values
              return(tcrossprod(y %*% diag(fn(v)), y))
            }

            ## ginvx is a helper to compute reciprocals
            ginvx <- function(x) {ifelse(x == 0,0,1/x)}

            ## decomposition for (a,b)
            gevd <- function(a, b = diag(nrow(a)), ncomps = 10) {
              bs <- mfunc(b,function(x) ginvx(sqrt(x)))
              #ev <- eigen(bs%*%a%*%bs)
              ev <- eigs(bs %*% a %*% bs, k = ncomps)
              return(list(gvalues = ev$values, gvectors = bs %*% ev$vectors))
            }

            x <- as.matrix(x)
            y <- as.matrix(y)

            if (!(nrow(x) == nrow(y)))
              stop("Number of rows in x, y matrixes is not equal")

            if (!is(kernel,"kernel"))
            {
              if (is(kernel,"function")) {
                kernel1 <- deparse(substitute(kernel))
                kernel2 <- deparse(substitute(kernel))
              }
              kernel1 <- do.call(kernel, kpar)
              kernel2 <- do.call(kernel, kpar2)
            }
            if (!is(kernel1,"kernel")) stop("kernel must inherit from class `kernel'")
            if (!is(kernel2,"kernel")) stop("kernel must inherit from class `kernel'")

            if (!is.null(w)) {
            w2 <- w/sum(w)
            cWps <- sqrt(w2)
            ## Generate LH
            Kx <- kernelMatrix(kernel1,x)
            Ky <- kernelMatrix(kernel2,y)
            n <- dim(Kx)[1]
            m <- 2
            VK <- matrix(0,n*2,n);
            VK[0:n,] <- t(cWps * t(Kx))
            VK[(n + 1):(2*n),] <- t(cWps*t(Ky))
            } else {
            Kx <- kernelMatrix(kernel1,x)
            Ky <- kernelMatrix(kernel2,y)
            n <- dim(Kx)[1]
            m <- 2
            ## Generate LH
            VK <- matrix(0,n*2,n);
            VK[0:n,] <- Kx
            VK[(n + 1):(2*n),] <- Ky
            }

            LH <- tcrossprod(VK, VK)

            for (i in 1:m)
              LH[((i - 1)*n + 1):(i*n),((i - 1)*n + 1):(i*n)] <- 0

            if (!is.null(w)) {
              RH <- matrix(0,n*m,n*m)
              RH[1:n,1:n] <- (t(t(Kx)*w2) + diag(rep(gamma,n))) %*% (Kx) + diag(rep(1e-6,n))
              RH[(n + 1):(2*n),(n + 1):(2*n)] <- (t(t(Ky)*w2) + diag(rep(gamma,n))) %*% (Ky) + diag(rep(1e-6,n))
              RH <- (RH + t(RH))/2
              ei <- gevd(LH, RH, ncomps = ncomps)
            } else {
              RH <- matrix(0,n*m,n*m)
              RH[1:n,1:n] <- (Kx + diag(rep(gamma,n))) %*% Kx + diag(rep(1e-6,n))
              RH[(n + 1):(2*n),(n + 1):(2*n)] <- (Ky + diag(rep(gamma,n))) %*% Ky + diag(rep(1e-6,n))
              RH <- (RH + t(RH))/2
              ei <- gevd(LH, RH, ncomps = ncomps)
            }
            ret <- new("kcca")

            ret@kcor <- as.double(ei$gvalues[1:ncomps])
            ret@xcoef <- matrix(as.double(ei$gvectors[1:n,1:ncomps]),n)
            ret@ycoef <- matrix(as.double(ei$gvectors[(n + 1):(2*n),1:ncomps]),n)
            ## xvar(ret) <- rotated(xpca) %*% cca$xcoef
            ## yvar(ret) <- rotated(ypca) %*% cca$ycoef
            return(ret)
          })

## gevd compute the generalized eigenvalue
## decomposition for (a,b)
.gevd <- function(a,b = diag(nrow(a)), ncomps=10) {
  bs <- .mfunc(b, function(x) .ginvx(sqrt(x)))
  #ev<-eigen(bs %*% a %*% bs)
  ev <- eigs(bs %*% a %*% bs, k = ncomps)
  return(list(gvalues = ev$values, gvectors = bs %*% ev$vectors))
}

## mfunc is a helper to compute matrix functions
.mfunc <- function(a,fn=sqrt) {
  e <- eigen(a); y <- e$vectors; v <- e$values
  return(tcrossprod(y %*% diag(fn(v)), y))
}

## ginvx is a helper to compute reciprocals
.ginvx <- function(x) {ifelse(x == 0, 0, 1/x)}
