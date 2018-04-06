#' A function to retrieve the BinaryID tables associated with a SpatialStreamNetwork
#' 
#'@description
#'
#'Retrieves the Binary ID tables which are stored in a database in the path slot of a SpatialStreamNetwork object and returns them as a list of data.frame objects.
#'
#'@usage
#'
#'\code{getBIDtables(ssn.object, networks = 1:nnetwork(ssn.object))}
#'
#'@param ssn.object An object of class SpatialStreamNetwork.
#'@param networks A numeric vector of the networks in a SpatialStreamNetwork object for which the binary ID tables should be extracted.
#'@return A named list of length \code{length(networks)} containing the binary ID tables for each of these networks. 
#'
#'@details
#'
#'The ssn.object must have a valid folder path stored in its path slot. This folder must contain a binaryID.db database.
#'
#'@export
getBIDtables <- function(
  ssn.object, 
  networks = 1:nnetwork(ssn.object)
  ){
  
  # Connect to database
  
  drvr <- dbDriver("SQLite")
  conn <- dbConnect(drvr, paste(ssn.object@path, "binaryID.db", sep = "/"))
  
  # Loop through database tables to get BID tables
  
  n.iter <- length(networks)
  tables <- vector("list", n.iter)
  
  for (i in 1:n.iter) {
    
    netid <- networks[i]
    binary_id_table <- dbReadTable(conn, paste0("net", netid))
    tables[[i]] <- binary_id_table[order(as.numeric(binary_id_table$rid)), ]
    
  }
  
  # disconnect from database
  
  dbDisconnect(conn)
  
  # return BID tables
  
  names(tables) <- paste0("net", networks)
  
  return(tables)
  
}
