% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{filterFeatures}
\alias{filterFeatures}
\title{Filter features}
\usage{
filterFeatures(data, feature, which.feature.type)
}
\arguments{
\item{data}{data frame}

\item{feature}{Character indicating which feature (i.e., data.frame column) to consider for filtering}

\item{which.feature.type}{Character vector indicating which feature types to include. If NULL, all are included.}
}
\value{
filtered data frame
}
\description{
This function filters a dataframe by a specified feature type(s)
}
