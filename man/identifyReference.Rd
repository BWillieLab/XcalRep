% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration_functions.R
\name{identifyReference}
\alias{identifyReference}
\title{Identify Reference Site}
\usage{
identifyReference(
  object,
  which.assay = NULL,
  which.parameter = NULL,
  which.phantom = NULL,
  which.data = "uncalibrated",
  pool.phantoms = F,
  return.scores = F
)
}
\arguments{
\item{object}{Calibration Object}

\item{which.assay}{Character specifying which assay to identify reference site for.}

\item{which.parameter}{Character specifying which parameters to consider when selecting reference site}

\item{which.phantom}{Character specifying which phantoms to select reference sites for. If unspecified, all phantoms are considered.}

\item{which.data}{Character specifying which data to identify reference site for.}

\item{pool.phantoms}{Logical specifying whether to pool ranking scores across all phantoms (true), or to rank phantoms separately (false, default).}

\item{return.scores}{Logical indicating whether list of scores are returned. If false (default), reference sites are returned only.}
}
\value{
If return.scores is false, named vector of character(s) specfying reference site is returned. In case of tie, multiple sites are returned. If return.scores is true, list of site rankings is returned.
}
\description{
Sites are ranked according to mean-squared error calculated between site-specific parameters and overall median parameter values (pooled across all timePoints and sites), and top ranking site is returned.
}
\seealso{
\code{\link{consistencyAnalysis}}, \code{\link{consistencyPlot}}
}
