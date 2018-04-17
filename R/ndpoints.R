#' Calculate the number of design points on a network
#'
#'@description
#'   
#'\code{ndpoints()} is a function that returns the number of design points on a SpatialStreamNetwork object, both in total and by network.
#'
#'@usage
#'
#'\code{ndpoints(ssn)}
#'
#'@param ssn An object of class SpatialStreamNetwork
#'@return A list with two elements: a vector of design points by network, and a scalar being the total number of design points in the SpatialStreamNetwork object.
#'
#'@details
#'
#'To find the number of prediction points on a SpatialStreamNetwork object, use \code{\link{nppoints}}.
#' 
#'@export
ndpoints <- function(ssn){
  
  ## Get vector of networks
  
  nets <- ssn@obspoints@SSNPoints[[1]]@point.data$netID
  nets <- as.numeric(as.character(nets))
  
  ## Get unique networks
  
  nets.u <- unique(nets)
  nets.u <- nets.u[order(nets.u)]
  n <- length(nets.u)
  
  ## Find number of points
  
  ndp <- vapply(nets, function(x){sum(nets == x)}, vector("numeric", n))
  names(ndp) <- paste("Net", nets.u)
  
  ## Return
  
  return(list(
    byNetwork = ndp,
    inTotal = sum(ndp)
  )
  )
  
}
