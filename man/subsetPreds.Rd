% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsetPreds.R
\name{subsetPreds}
\alias{subsetPreds}
\title{Subset the prediction sites in a SpatialStreamNetwork}
\usage{
subsetPreds(ssn, subset)
}
\arguments{
\item{ssn}{An object of class SpatialStreamNetwork}

\item{subset}{A logical condition used to subset the prediction sites.}
}
\value{
A SpatialStreamNetwork. The new preds shapefile overwrites the old one in \code{ssn@path}.
}
\description{
This function works like subsetSSN, but instead of focussing on the observed sites, this allows the prediction sites alone to be subsetted according to some logical criterion or criteria. It is particularly useful when attempting to keep only prediction sites with certain PIDs.
}
