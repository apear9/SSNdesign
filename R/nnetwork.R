#' Find the number of distinct networks in a SpatialStreamNetwork
#'
#'@description
#'
#'Calculates the number of networks are present inside of a SpatialStreamNetwork object. 
#'
#'@param ssn an object of class SpatialStreamNetwork
#'@return a numeric scalar
#'
#'@examples
#'
#'s <- createSSN(c(10, 10), binomialDesign(c(2, 2)), path = paste(tempdir(), "tmp.ssn", sep = "/"), importToR = TRUE)
#'nnetwork(s) # reports 2
#'      
#'@export
nnetwork <- function(ssn){
  if(class(ssn)[1] != "SpatialStreamNetwork"){
    stop("Please provide an object of class SpatialStreamNetwork")
  }
  networks <- unique(ssn@network.line.coords$NetworkID)
  nnetworks <- length(networks)
  return(nnetworks)
}