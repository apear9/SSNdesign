#' Calculate Shreve stream order and afv
#' 
#'@description 
#'
#' Calculates Shreve stream order and additive function values on an arbitrary SpatialStreamNetwork object.
#' 
#'@usage 
#' 
#'\code{calculateShreveStreamOrderAndAFVs(ssn)}
#' 
#'@param ssn an object of class SpatialStreamNetwork
#'@return An object of class SpatialStreamNetwork. The SSNPoints for the obspoints slot will be updated to reflect the selected design. 
#'
#'@details 
#'
#'This function accepts a SpatialStreamNetwork object and calculates Shreve stream orders and additive function values for it, which are then joined back to the `ssn@data` slot. This funciton should be used when a SpatialStreamNetwork has no edge attributes which could be used to calculate the AFV. 
#'    
#'@export
#'
calculateShreveStreamOrderAndAFVs <- function(ssn, BID.tables = NULL){
  
  # Calculate how many networks there are
  nn <- nnetwork(ssn)
  
  # Basic checking of BID.tables input
  if(!is.null(BID.tables)){
    if(!is.list(BID.tables) | !(length(BID.tables) == nn)){
      stop("Please ensure that BID.tables is a list of the binary ID tables for all networks in this SpatialStreamNetwork object.")
    }
  }
  
  # How many segments
  ns <- length(ssn@data$rid)
  
  # Extract out the binary ID tables for this spatial stream network object
  if(is.null(BID.tables)){
    BID.tables <- getBIDtables(ssn)
    if(ns > 1e4){
      print("Large spatial stream network detected. Reading in the binary ID tables may take several minutes.")
    }
  }
  print("Hi")
  # Calculate Shreve stream order per network
  shreve.orders <- vector("list", nn)
  for(i in 1:nn){
    shreve.orders[[i]] <- calculateShreveStreamOrder(BID.tables[[i]])
  }
  print("Hi again")
  # Calculate additive function values
  afv.networks <- lapply(
    shreve.orders, 
    function(x){
      afv <- x$shreve/max(x$shreve)
      return(data.frame(
        rid = x$rid, 
        shreve = x$shreve,
        afv = afv
      )
      )
    }
  )
  print("Hi lmao")
  # Join each to a part of the ssn@data table
  ssn@data$shreve <- 0
  ssn@data$afv <- 0
  for(i in 1:nn){
    tmp <- afv.networks[[i]]
    ind <- match(tmp$rid, ssn@data$rid)
    ssn@data$shreve[ind] <- tmp$shreve
    ssn@data$afv[ind] <- tmp$afv
  }

  # Return output
  return(ssn)
  
}
