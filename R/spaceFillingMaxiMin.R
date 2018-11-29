#' A utility function for maximin space-filling designs
#' 
#'@description
#'   
#'\code{spaceFillingMaximin} is a utility function used by \code{\link{constructSpaceFillingDesign}}. It is also fully compatible with \code{\link{optimiseSSNDesign}}. It is a utility function that maximises the minimum interpoint distance between a set of design points.  
#' 
#'@usage
#'
#'\code{spaceFillingMaximin(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
#' 
#'@param ssn An object of class SpatialStreamNetwork. This argument may be ignored, but it is present so this function is consistent with all other utility functions.
#'@param glmssn An model object of class glmssn. Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#'@param prior.parameters A list of random functions that are parameterised in terms of n.draws.
#'@param n.draws A numeric scalar for the number of Monte Carlo draws to use when approximating the utility. Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{optimiseSSNDesign}}.
#'@return A numeric scalar.
#' 
#'@details
#'
#'\code{spaceFillingMaximin} is deterministic and will ignore the arguments \code{prior.parameters} and \code{n.draws}. It is expected to be very fast, unless there are several hundred design points to select. 
#' 
#'@export
spaceFillingMaxiMin <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  d <- extra.arguments$Matrices.Obs$d
  ij <- row.names(d) %in% design.points
  d <- d[ij, ij]
  v <- d[upper.tri(d)] 
  
  return(min(v[v != 0])[1]) 
  
}
