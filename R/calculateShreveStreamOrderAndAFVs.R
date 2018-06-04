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
#'@param BID.tables optionally provide pre-computed binary ID tables. This can be beneficial if the stream network is large or if the binary ID tables have already been computed for another purpose.
#'@return An object of class SpatialStreamNetwork. 
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

  # Calculate Shreve stream order per network
  shreve.orders <- vector("list", nn)
  for(i in 1:nn){
    shreve.orders[[i]] <- calculateShreveStreamOrder(BID.tables[[i]])
  }
 
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

  # Join each to a part of the ssn@data table
  ssn@data$shreve <- 0
  ssn@data$afv <- 0
  for(i in 1:nn){
    tmp <- afv.networks[[i]]
    ind <- match(tmp$rid, ssn@data$rid)
    ssn@data$shreve[ind] <- tmp$shreve
    ssn@data$afv[ind] <- tmp$afv
  }
  
  # Extract shreve and additive function value columns to the obspoints and predpoints frames
  ssn <- extractShreveAndAFVs(ssn)

  # Return output
  return(ssn)
  
}
