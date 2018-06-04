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
#'@return An object of class SpatialStreamNetwork. 
#'
#'@details
#'
#'Using this function directly is not recommended.
#'
#'@export
extractShreveAndAFVs <- function(ssn){
  
  # Check whether prediction sites exist
  preds.exist <- length(ssn@predpoints@SSNPoints) > 0
  # Some important numbers
  # # Add shreve and afv columns to observed +/- pred points and join
  ssn@obspoints@SSNPoints[[1]]@point.data$shreve <- 0
  ssn@obspoints@SSNPoints[[1]]@point.data$afv <- 0
  obs.unique.rids <- unique(ssn@obspoints@SSNPoints[[1]]@point.data$rid)
  obs.unique.rids.ind <- match(obs.unique.rids, ssn@data$rid)
  obs.unique.shreve <- ssn@data$shreve[obs.unique.rids.ind]
  obs.unique.afv <- ssn@data$afv[obs.unique.rids.ind]
  n.unique.rids <- length(obs.unique.rids)
  for(i in 1:n.unique.rids){
    ind <- which(ssn@obspoints@SSNPoints[[1]]@point.data$rid == obs.unique.rids[i])
    ssn@obspoints@SSNPoints[[1]]@point.data$shreve[ind] <- obs.unique.shreve[i]
    ssn@obspoints@SSNPoints[[1]]@point.data$afv[ind] <- obs.unique.afv[i]
  }
  if(preds.exist){
    ssn@predpoints@SSNPoints[[1]]@point.data$shreve <- 0
    ssn@predpoints@SSNPoints[[1]]@point.data$afv <- 0
    pred.unique.rids <- unique(ssn@predpoints@SSNPoints[[1]]@point.data$rid)
    pred.unique.rids.ind <- match(pred.unique.rids, ssn@data$rid)
    pred.unique.shreve <- ssn@data$shreve[pred.unique.rids.ind]
    pred.unique.afv <- ssn@data$afv[pred.unique.rids.ind]
    n.unique.rids <- length(pred.unique.rids)
    for(i in 1:n.unique.rids){
      ind <- which(ssn@predpoints@SSNPoints[[1]]@point.data$rid == pred.unique.rids[i])
      ssn@predpoints@SSNPoints[[1]]@point.data$shreve[ind] <- pred.unique.shreve[i]
      ssn@predpoints@SSNPoints[[1]]@point.data$afv[ind] <- pred.unique.afv[i]
    }
  }
  
  # Return
  return(ssn)
  
}