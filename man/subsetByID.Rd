% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsetSSNbyID.R
\name{subsetByID}
\alias{subsetByID}
\title{Subset the observed and/or prediction sites by ID}
\usage{
subsetByID(ssn, new.ssn.path, obs, preds, locID = FALSE)
}
\arguments{
\item{ssn}{An object of class SpatialStreamNetwork.}

\item{new.ssn.path}{A path to a new .ssn directory where the outputs of this function should be stored.}

\item{obs}{A numeric or character vector containing the pids or locIDs of the observed sites that should be kept in the subset.}

\item{preds}{A numeric or character vector containing the pids or locIDs of the prediction sites that should be kept in the subset. This argument can be skipped.}

\item{locID}{A logical indicating whether the obs and preds vectors contain pids or locIDs. Defaults to FALSE.}
}
\value{
An object of class SpatialStreamNetwork containing only the specified subset of sites.
}
\description{
This is a function to subset the observed and/or prediction sites on a SpatialStreamNetwork based on pids or locIDs.
}
\details{
Note, the SpatialStreamNetwork that is returned will have been re-imported from the \code{new.ssn.path} directory after the subset operations are complete.
}
