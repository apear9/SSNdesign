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
