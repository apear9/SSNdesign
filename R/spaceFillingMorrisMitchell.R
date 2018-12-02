#' A utility function for maximin space-filling designs per the criterion proposed by Morris and Mitchell (1995)
#' 
#'@description
#'   
#'\code{spaceFillingMorrisMitchell} is a utility function used by the function \code{\link{constructSpaceFillingDesign}}. It can also be used with \code{\link{optimiseSSNDesign}}. It is a utility function that maximises the minimum interpoint distance between a set of design points.  
#' 
#'@param ssn Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param glmssn Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#'@param prior.parameters Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param n.draws Ignored in this function but this argument is present so this function is consistent with all other utility functions.
#'@param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{optimiseSSNDesign}}.
#'@return The expected utility of a design under the Morris and Mitchell (1995) version of the maximin space-filling utility function.
#' 
#'@details
#'
#'\code{spaceFillingMorrisMitchell} is deterministic and will ignore the arguments \code{prior.parameters} and \code{n.draws}. It is expected to be very fast, unless there are several hundred design points to select. 
#'
#'For examples, please see \code{\link{constructSpaceFillingDesign}}.
#' 
#' @references 
#' 
#' Morris, M.D. & Mitchell, T.J. (1995). Exploratory Designs for Computational Experiments. \emph{Journal of Statistical Planning and Inference}, \emph{43}, 381-402. 
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
