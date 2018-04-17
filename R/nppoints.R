#' Calculate the number of prediction points on a network
#'
#'@description
#'   
#'\code{nppoints()} is a function that returns the number of prediction points on a SpatialStreamNetwork object, both in total and by network.
#'
#'@usage
#'
#'\code{nppoints(ssn)}
#'
#'@param ssn An object of class SpatialStreamNetwork
#'@return A list with two elements: a vector of prediction points by network, and a scalar being the total number of prediction points in the SpatialStreamNetwork object.
#'
#'@details
#'
#'To find the number of design points on a SpatialStreamNetwork object, use \code{\link{ndpoints}}.
#' 
#'@export
nppoints <- function(ssn){
  
  ## Get vector of networks
  
  nets <- ssn@predpoints@SSNPoints[[1]]@point.data$netID
  nets <- as.numeric(as.character(nets))
  
  ## Get unique networks
  
  nets.u <- unique(nets)
  nets.u <- nets.u[order(nets.u)]
  
  ## Find number of points
  
  npp <- vapply(nets.u, function(x){sum(nets == x)}, vector("numeric", 1))
  names(ndp) <- paste("Net", nets.u)
  
  ## Return
  
  return(list(
    byNetwork = npp,
    inTotal = sum(npp)
  )
  )
  
}
