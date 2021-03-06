% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixing_utilities.R
\name{m3_mixing_model}
\alias{m3_mixing_model}
\title{m3_mixing_model}
\usage{
m3_mixing_model(gdata, end.members = NULL)
}
\arguments{
\item{gdata}{A geochemical_dataset object with the major ions only.}

\item{end.members}{If provided, this is a numeric vector with the indices of the end members
to be used in the mixing calculations. If it is not specified then the end members are
defined using the convex hull of the first two principal components.}
}
\value{
This function returns a list with the following entries:
\itemize{
\item mixing.ratios: Matrix with the mixing ratios estimated using constrained
least-squares.
\item res.pca: Object returned by the princomp function. This object contains the
eigenvectors, eigenvalues, scores and other results of the PCA.
\item res.lm: List with the results of the multiple linear regression applied on each
ion included in the gdata geochemical_dataset object. Each entry of this list is a lm
object for each ion.
\item names.lm: Character vector with the names of the ions used in the mixing analysis.
\item X:
\item dataset:
}
}
\description{
Function that implements the M3 mixing model. This model is based on the projection of
the chemical information in a low-dimensional space defined using Principal Component
Analysis (PCA). The compositions of pre-defined end-members are also projected and from
these coordinates, the mixing ratios of each end-member are calculated using a constrained
least-squares procedure (mixing ratios add up to 1).
}
\references{
Laaksoharju, M., Skarrman, C., \& Skarrman, E. (1999). Multivariate mixing and mass
balance (M3) calculations, a new tool for decoding hydrogeochemical information.
Applied Geochemistry, 14(7), 861–871. http://doi.org/10.1016/S0883-2927(99)00024-4
}
\seealso{
Other mixing functions: \code{\link{constrained_lm}},
  \code{\link{mix_model}},
  \code{\link{plot_m3_mixing_results}},
  \code{\link{select_end_members}},
  \code{\link{uncertainty_m3_mixing_model}}
}
\author{
Oscar Garcia-Cabrejo \email{khaors@gmail.com}
}
\concept{mixing functions}
