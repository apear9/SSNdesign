#' Extract Shreve stream order and additive function values to observed and prediction points on a SpatialStreamNetwork
#' 
#'@description 
#'
#' This function is designed to be used internally to join the appropriate shreve stream order and additive function values to the observed and prediction points on a SpatialStreamNetwork when \code{calculateShreveStreamOrderAndAFVs()} is run.
#' 
#'@usage 
#' 
#'\code{extractShreveAndAFVs(ssn)}
#' 
#'@param ssn an object of class SpatialStreamNetwork
#'@param afv.networks a list that is generated inside \code{calculateShreStreamOrderAndAFVs()}
#'@return An object of class SpatialStreamNetwork. 
#'
#'@details
#'
#'Using this function directly is not recommended.
#'
#'@export
extractShreveAndAFVs(ssn, afv.networks){
  
  # Check whether prediction sites exist
  preds.exist <- length(ssn@predpoints@SSNPoints) > 0
  # Check whether edges have Shreve and afv columns
  if(is.null(ssn@data$shreve) | is.null(ssn@data$afv)){
    stop("This function cannot work unless columns called 'shreve' and 'afv' exist inside ssn@data")
  }
  # Some important numbers
  nn <- length(afv.networks) # num of networks
  # Add shreve and afv columns to observed +/- pred points 
  ssn@obspoints@SSNPoints[[1]]@point.data$shreve <- 0
  ssn@obspoints@SSNPoints[[1]]@point.data$afv <- 0
  if(preds.exist){
    ssn@predpoints@SSNPoints[[1]]@point.data$shreve <- 0
    ssn@predpoints@SSNPoints[[1]]@point.data$afv <- 0
  }
  # Perform the extraction
  for(i in 1:nn){
    # Observed
    any.in <- sum(ssn@obspoints@SSNPoints[[1]]@point.data$rid %in% afv.networks[[i]]$rid) != 0
    if(any.in){
      ind.obs.to <- match(ssn@obspoints@SSNPoints[[1]]@point.data$rid, afv.networks[[i]]$rid)
      ind.obs.from <- which(ssn@obspoints@SSNPoints[[1]]@point.data$rid %in% afv.networks[[i]]$rid)
      ssn@obspoints@SSNPoints[[1]]@point.data$shreve[ind.obs.from] <- afv.networks[[i]]$shreve[ind.obs.to]
      ssn@obspoints@SSNPoints[[1]]@point.data$afv[ind.obs.from] <- afv.networks[[i]]$afv[ind.obs.to]
    }
    # To prediction, if relevant
    if(preds.exist){
      any.in <- sum(ssn@predpoints@SSNPoints[[1]]@point.data$rid %in% afv.networks[[i]]$rid) != 0
      if(any.in){
        ind.prd.to <- match(ssn@predpoints@SSNPoints[[1]]@point.data$rid, afv.networks[[i]]$rid)
        ind.prd.from <- which(ssn@predpoints@SSNPoints[[1]]@point.data$rid %in% afv.networks[[i]]$rid)
        ssn@predpoints@SSNPoints[[1]]@point.data$shreve[ind.prd.from] <- afv.networks[[i]]$shreve[ind.prd.to]
        ssn@predpoints@SSNPoints[[1]]@point.data$afv[ind.prd.from] <- afv.networks[[i]]$afv[ind.prd.to]
      }
    }
  }
  
  # Return
  return(ssn)
  
}