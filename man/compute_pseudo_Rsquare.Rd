% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.compute_pseudo_Rsquare.R
\name{compute_pseudo_Rsquare}
\alias{compute_pseudo_Rsquare}
\title{ccompute_pseudo_Rsquare}
\usage{
compute_pseudo_Rsquare(
  ACTUAL,
  PRED,
  wt,
  type = c("mean", "sum"),
  ACTUAL_mean = FALSE
)
}
\arguments{
\item{ACTUAL}{a numeric matrix, actual observations.}

\item{PRED}{a numeric matrix, predictions.}

\item{wt}{a numeric vector, sample weights.}

\item{type}{a string, type of error measure, (either \code{mean} or \code{sum}).}

\item{ACTUAL_mean}{if TRUE use  the actual mean in the numerator if FALSE use the predicted mean}
}
\value{
A double, the value of the objective function.
}
\description{
Compute the value of the compute_pseudo_Rsquare ans make the mean or the sum
}
\details{
The standardized version of the objective function, computed using \code{type = "wsRMSE"}, enables us to deal with the fact the continuous variables are measured on different scales.
}
