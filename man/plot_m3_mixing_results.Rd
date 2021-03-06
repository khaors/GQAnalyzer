% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixing_utilities.R
\name{plot_m3_mixing_results}
\alias{plot_m3_mixing_results}
\title{plot_m3_mixing_results}
\usage{
plot_m3_mixing_results(gdata, mixing.res, type = c("concentration",
  "mixing.ratio", "residual"), element)
}
\arguments{
\item{gdata}{A geochemical_dataset object}

\item{mixing.res}{A list with ther results of the m3_mixing_model function}

\item{type}{A character string specifying the plot type to be created. Currently the
supported options are:
\itemize{
\item concentration
\item mixing.ratio
}}

\item{element}{A character string with the name of the ion}
}
\value{
This function returns a ggplot2 object with the requested plot.
}
\description{
Function to plot different results obtained from the application of the M3 mixing model
}
\seealso{
Other mixing functions: \code{\link{constrained_lm}},
  \code{\link{m3_mixing_model}}, \code{\link{mix_model}},
  \code{\link{select_end_members}},
  \code{\link{uncertainty_m3_mixing_model}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
\concept{mixing functions}
