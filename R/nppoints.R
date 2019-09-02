#' Calculate the number of prediction points on a network
#'
#'@description
#'   
#' The function \code{nppoints} returns the number of prediction points on a SpatialStreamNetwork object, both in total and in each network.
#'
#'@param ssn An object of class SpatialStreamNetwork
#'@param use A string indicating whether the pids or locIDs should be counted. 
#'@return A list with two elements: a vector of prediction points by network, and a scalar being the total number of prediction points in the SpatialStreamNetwork object.
#'
#'@details
#'
#'To find the number of design points on a SpatialStreamNetwork object, use \code{\link{ndpoints}}.
#'
#'@examples
#'
#'set.seed(1)
#'
#'s1 <- createSSN(10, binomialDesign(10),binomialDesign(10), 
#'path = paste(tempdir(), "s1_no_reps_preds.ssn", sep = "/"), importToR = TRUE)
#'nppoints(s1)
#'
#'s2 <- createSSN(10, binomialDesign(10, 2, "Time"), binomialDesign(10, 2, "Time"),
#' path = paste(tempdir(), "s2_reps_preds.ssn", sep = "/"), importToR = TRUE)  
#'nppoints(s2) # total number of observations
#'nppoints(s2, "locIDs") # total number of sites
#' 
#'@export
nppoints <- function(ssn, use = "pids"){
  
  ## Check whether prediction points exist
  preds.exist <- length(ssn@predpoints@SSNPoints) > 0
  if(!preds.exist){
    warning("There are no prediction points in this SpatialStreamNetwork object.")
    return(0)
  }
  
  # Early exits for bad inputs
  if(!(use %in% c("pids", "locIDs"))){
    stop("The argument use must be either pids or locIDs. Check spelling.")
  }
  
  if(use == "pids"){
    ## Get vector of networks
    
    nets <- ssn@predpoints@SSNPoints[[1]]@point.data$netID
    nets <- as.numeric(as.character(nets))
    
    ## Get unique networks
    
    nets.u <- unique(nets)
    nets.u <- nets.u[order(nets.u)]
    
    ## Find number of points
    
    # Per network
    npp.net <- vapply(nets.u, function(x){sum(nets == x)}, vector("numeric", 1))
    names(npp.net) <- paste("Net", nets.u)
    
    # Overall
    npp <- sum(npp.net)
    
  } else {
    
    ## Get information on unique locations
    
    locs.u <- unique(ssn@predpoints@SSNPoints[[1]]@point.data$locID)
    
    ## Per network
    
    locs.nets.u <- unique(ssn@predpoints@SSNPoints[[1]]@point.data[,c("locID", "netID")])
    nets.u <- unique(locs.nets.u$netID)
    nets.n <- length(nets.u)
    
    ## Summarise information
    
    # Per network
    
    npp.net <- vector("numeric", nets.n)
    for(i in 1:nets.n){
      npp.net[i] <- sum(locs.nets.u$netID == nets.u[i])
    }
    names(npp.net) <- paste("Net", nets.u)
    
    # Overall
    npp <- length(locs.u)
    
  }
  
  ## Return
  
  return(list(
    byNetwork = npp.net,
    inTotal = npp
  )
  )
  
}
