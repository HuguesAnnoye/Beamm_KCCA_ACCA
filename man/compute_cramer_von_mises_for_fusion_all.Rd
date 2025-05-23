% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.compute_cramer_von_mises_for_fusion_all.R
\name{compute_cramer_von_mises_for_fusion_all}
\alias{compute_cramer_von_mises_for_fusion_all}
\title{cramer von minsen discrete}
\usage{
compute_cramer_von_mises_for_fusion_all(
  ACTUAL,
  PRED,
  wt = NULL,
  wt_PRED = NULL,
  names.NCV = NULL,
  type = c("mean", "sum", "sumst"),
  parallel = FALSE,
  nnodes = 10
)
}
\arguments{
\item{ACTUAL}{a data frame, actual observations.}

\item{PRED}{a data frame, predictions.}

\item{wt}{a numeric vector, sample weights for ACTUAL.}

\item{wt_PRED}{a numeric vector, sample weights for ACTUAL.}

\item{names.NCV}{a character vector with the names of the non-common variables used for the prediction.}

\item{type}{a string, type of error measure, (either \code{mean}, \code{sum} or \code{sumst}).}

\item{parallel}{a logical, whether or not use parallel computing.}

\item{nnodes}{an integer, number of parallel session (CPUs).}
}
\value{
A list,with the values of the objective function.
}
\description{
Compute the value of the compute_pseudo_Rsquare ans make the mean or the sum
}
