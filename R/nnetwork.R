#' Find the number of distinct networks in a SpatialStreamNetwork
#' 
#' \code{nnetwork()} is a function to find how many separate networks there are in a single SpatialStreamNetwork object.
#' 
#' @param ssn an object of class SpatialStreamNetwork
#' @return a single number
#' 
#' @examples
#' \dontrun{NONE AS YET}
#' 
#' @export
nnetwork <- function(ssn){
  if(class(ssn)[1] != "SpatialStreamNetwork"){
    stop("Please provide an object of Class SpatialStreamNetwork")
  }
  networks <- unique(ssn@network.line.coords$NetworkID)
  nnetworks <- length(networks)
  return(nnetworks)
}