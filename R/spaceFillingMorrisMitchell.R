#' A utility function for maximin space-filling designs per the criterion proposed by Morris and Mitchell (1995)
#' 
#'@description
#'   
#'\code{spaceFillingMorrisMitchell} is a utility function that can be used with \code{\link{findOptimalDesign}}. It is a utility function that maximises the minimum interpoint distance between a set of design points.  
#' 
#'@usage
#'
#'\code{spaceFillingMorrisMitchell(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
#' 
#'@param ssn An object of class SpatialStreamNetwork. This argument may be ignored, but it is present so this function is consistent with all other utility functions.
#'@param glmssn An model object of class glmssn. Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#'@param prior.parameters A list of random functions that are parameterised in terms of n.draws.
#'@param n.draws A numeric scalar for the number of Monte Carlo draws to use when approximating the utility. Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{findOptimalDesign}}.
#'@return A numeric scalar.
#' 
#'@details
#'
#'\code{spaceFillingMorrisMitchell} is deterministic and will ignore the arguments \code{prior.parameters} and \code{n.draws}. It is expected to be very fast, unless there are several hundred design points to select. 
#' 
#'@export
spaceFillingMorrisMitchell <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  d <- extra.arguments$Matrices.Obs$d
  p <- extra.arguments$p
  ij <- row.names(d) %in% design.points
  d <- d[ij, ij]
  v <- d[upper.tri(d) & d != 0] 
  ds <- sort(unique(v))
  dsp <- ds^-p
  J <- numeric(length(ds))
  for(i in 1:length(ds)){
    J[i] <- sum(ds == ds[i]) 
  }
  
  return(-sum((J * dsp))^(1/p)) 
  
}
