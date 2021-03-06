% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractStreamEdgeCovariates.R
\name{extractStreamEdgeCovariates}
\alias{extractStreamEdgeCovariates}
\title{Extract covariates on stream edges to observed and/or predicted sites}
\usage{
extractStreamEdgeCovariates(ssn, columns)
}
\arguments{
\item{ssn}{an object of class SpatialStreamNetwork}

\item{columns}{a vector of column names that should be extracted from the data slot of the SpatialStreamNetwork object}
}
\value{
An object of class SpatialStreamNetwork.
}
\description{
This function transfers specified variables recorded on the edges of a SpatialStreamNetwork object to the observed (and possibly prediction) sites in that SpatialStreamNetwork object.
}
\examples{

set.seed(1)

# Create an SSN
s <- createSSN(10, binomialDesign(10), 
path = paste(tempdir(), "example01.ssn", sep = "/"), importToR = TRUE)

# Extract stream edge covariates
s <- extractStreamEdgeCovariates(s, "Length")

}
