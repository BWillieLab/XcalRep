% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_functions.R
\name{analyzeWhich}
\alias{analyzeWhich}
\title{Specify features to analyze}
\usage{
analyzeWhich(
  object,
  include.sites = "all",
  include.sections = "all",
  include.times = "all",
  include.phantoms = "all",
  include.parameters = "all",
  include.scanDates = "all",
  include.scanIDs = "all",
  omit.sites = NULL,
  omit.sections = NULL,
  omit.times = NULL,
  omit.phantoms = NULL,
  omit.parameters = NULL,
  omit.scanDates = NULL,
  omit.scanIDs = NULL,
  which.assay = NULL
)
}
\arguments{
\item{object}{calibration object}

\item{include.sites}{sites to include in analysis}

\item{include.sections}{sections to include in analysis}

\item{include.times}{timePoints to include in analysis}

\item{include.phantoms}{phantoms to include in analysis}

\item{include.parameters}{parameters to include in analysis}

\item{include.scanDates}{scanDates to include in analysis}

\item{include.scanIDs}{scanIDs to include in analysis}

\item{omit.sites}{sites to omit from analysis}

\item{omit.sections}{phantom sections to omit from analysis}

\item{omit.times}{timePoints to omit from analysis}

\item{omit.phantoms}{phantoms to omit from analysis}

\item{omit.parameters}{parameters to omit from analysis}

\item{omit.scanDates}{scanDates to omit from analysis}

\item{omit.scanIDs}{scanIDs to omit from analysis}

\item{which.assay}{specifies assay}
}
\value{
list
}
\description{
Preprocess HR-pQCT data for subsequent analysis.
Inclusion parameters are considered first, followed by omission parameters.
For inclusion parameters, "all", c("parameter_names") or NULL are supported inputs.
For omission parameters, c("parameter_names") or NULL are supported inputs.
}
