% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{calibrationPlot}
\alias{calibrationPlot}
\title{Plot Calibration Curves}
\usage{
calibrationPlot(
  object,
  which.assay = NULL,
  which.parameter = NULL,
  which.time = NULL,
  return.plt.handle = F,
  color.option = "cividis",
  color.begin = 0,
  color.end = 0,
  point.size = 2
)
}
\arguments{
\item{object}{Calibration Object.}

\item{which.assay}{Character specifying which assay to get calibration curves from.}

\item{which.parameter}{Vector specifying which parameter calibration curves to plot. If unspecified, all parameters plotted.}

\item{which.time}{Vector specifying which times to plot calibration curves for. If unspecified, all timePoints plotted.}

\item{return.plt.handle}{Logical specifying whether to return list of plot handles. If TRUE, list of plot handles is returned. If FALSE, plots are printed without returning handle.}

\item{color.option}{Character indicating which color map option to use (from viridis palette). One of:
\itemize{
\item "magma" (or "A")
\item "inferno" (or "B")
\item "plasma" (or "C")
\item "viridis" (or "D")
\item "cividis" (or "E") - Default
}}

\item{color.begin}{Hue in [0,1] at which the viridis color map begins. Default is 0.}

\item{color.end}{Hue in [0,1] at which the viridis color map ends. Default is 0.}

\item{point.size}{Numeric specifying size of scatter plot points. Only implemented if which.plot is "scatter". Default is 2.}
}
\value{
plot handle (see return.plt.handle argument)
}
\description{
Pair-wise calibration curves for each parameter at each timePoint are plotted with 95% confidence interval and dashed line reference indicating x=y.
}
\details{
Currently intercepts and slopes are fitted for all calibration curves, independent of whether intercepts were fit for calibration equation using fit.equation. This will be updated to accomodate no intercept calibration curve plots in future releasaed of XcalRep.
}
\seealso{
\code{\link{fitCalibration}}
}
