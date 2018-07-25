#' Calculate the number of design points on a network
#'
#'@description
#'   
#'\code{ndpoints()} is a function that returns the number of design points on a SpatialStreamNetwork object, both in total and by network.
#'
#'@usage
#'
#'\code{ndpoints(ssn, use = "pids")}
#'
#'@param ssn An object of class SpatialStreamNetwork
#'@param use A string indicating whether pids or locIDs should be used to count sites.
#'@return A list with two elements: a vector of design points by network, and a scalar being the total number of design points in the SpatialStreamNetwork object.
#'
#'@details
#'
#'To find the number of prediction points on a SpatialStreamNetwork object, use \code{\link{nppoints}}.
#' 
#'@export
ndpoints <- function(ssn, use = "pids"){
  
  # Early exits for bad inputs
  if(!(use %in% c("pids", "locIDs"))){
    stop("The argument use must be either pids or locIDs. Check spelling.")
  }
  
  if(use == "pids"){
    ## Get vector of networks
    
    nets <- ssn@obspoints@SSNPoints[[1]]@point.data$netID
    nets <- as.numeric(as.character(nets))
    
    ## Get unique networks
    
    nets.u <- unique(nets)
    nets.u <- nets.u[order(nets.u)]
    
    ## Find number of points
    
    # Per network
    ndp.net <- vapply(nets.u, function(x){sum(nets == x)}, vector("numeric", 1))
    names(ndp.net) <- paste("Net", nets.u)
    
    # Overall
    ndp <- sum(ndp.net)
    
  } else {
    
    ## Get information on unique locations
    
    locs.u <- unique(ssn@obspoints@SSNPoints[[1]]@point.data$locID)
    
    ## Per network
    
    locs.nets.u <- unique(ssn@obspoints@SSNPoints[[1]]@point.data[,c("locID", "netID")])
    nets.u <- unique(locs.nets.u$netID)
    nets.n <- length(nets.u)
    
    ## Summarise information
    
    # Per network
    
    ndp.net <- vector("numeric", nets.n)
    for(i in 1:nets.n){
      ndp.net[i] <- sum(locs.nets.u$netID == nets.u[i])
    }
    names(ndp.net) <- paste("Net", nets.u)
    
    # Overall
    ndp <- length(locs.u)
    
  }
  
  ## Return
  
  return(list(
    byNetwork = ndp.net,
    inTotal = ndp
  )
  )
  
}
