#' Find the number of distinct networks in a SpatialStreamNetwork
#'
#'@description
#'
#'Calculates the number of networks are present inside of a SpatialStreamNetwork object. 
#'
#'@usage
#'   
#'\code{nnetwork(ssn)}
#'
#'@param ssn an object of class SpatialStreamNetwork
#'@return a numeric scalar
#' 
#'@export
nnetwork <- function(ssn){
  if(class(ssn)[1] != "SpatialStreamNetwork"){
    stop("Please provide an object of Class SpatialStreamNetwork")
  }
  networks <- unique(ssn@network.line.coords$NetworkID)
  nnetworks <- length(networks)
  return(nnetworks)
}