nnetwork <- function(ssn){
  if(class(ssn)[1] != "SpatialStreamNetwork"){
    stop("Please provide an object of Class SpatialStreamNetwork")
  }
  networks <- unique(ssn@network.line.coords$NetworkID)
  nnetworks <- length(networks)
  return(nnetworks)
}