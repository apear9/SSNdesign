% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCalculatedSE.R
\name{getCalculatedSE}
\alias{getCalculatedSE}
\title{A convenience function that extracts the standard errors for the covariance parameters from a \code{glmssn} object}
\usage{
getCalculatedSE(glmssn, log.scale = FALSE)
}
\arguments{
\item{glmssn}{A fitted \code{glmssn} object.}

\item{log.scale}{Whether the standard errors should be on the log-scale. Defaults to \code{FALSE}.}
}
\value{
A vector of expected standard errors for the glmssn object's covariance parameters.
}
\description{
This function returns the expected standard errors on the covariance parameters for a given set of parameter values. It is a more limited version of \code{\link{getCalculatedCovMatrix}}. This function returns only the square roots of the diagonals of the observed fisher information matrix.
}
\details{
The transformation of the covariance matrix to and from the log-scale is performed using the delta method.
}
