#' A utility function for maximin space-filling designs
#' 
#'@description
#'   
#'\code{spaceFillingMaximin} is a utility function used by \code{\link{constructSpaceFillingDesign}}. It is also fully compatible with \code{\link{optimiseSSNDesign}}. It is a utility function that maximises the minimum interpoint distance between a set of design points.  
#' 
#'@param ssn Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param glmssn Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#'@param prior.parameters Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param n.draws Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{optimiseSSNDesign}}.
#'@return The expected utility of a design. This number represents the minimum (non-zero) distance between any pair of points in the design.
#' 
#'@details
#'
#'\code{spaceFillingMaximin} is deterministic and will ignore the arguments \code{prior.parameters} and \code{n.draws}. It is expected to be very fast, unless there are several hundred design points to select. 
#' 
#' For examples, please see \code{\link{constructSpaceFillingDesign}}.
#' 
#'@export
spaceFillingMaxiMin <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  d <- extra.arguments$Matrices.Obs$d
  ij <- row.names(d) %in% design.points
  d <- d[ij, ij]
  v <- d[upper.tri(d)] 
  
  return(min(v[v != 0])[1]) 
  
}
