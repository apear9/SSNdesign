writeStreams <- function(x, new.ssn.path, ...){
  
  # Check inputs 
  if(!isSSN(x)){
    stop("The argument x must be an object of class SpatialStreamNetwork.")
  }
  if(!dir.exists(new.ssn.path)){
    tryCatch(dir.create(new.ssn.path), error = function(e) stop("Check path specification. Could not create directory."))
  }
  
  # Copy database and netID files to new folder
  db <- dir(x@path, "binaryID", full.names = TRUE)
  db.success <- file.copy(db, new.ssn.path)
  if(!all(db.success)){
    warning("Database file may not have copied correctly. Try specifying a new ssn path and try again.")
  }
  
  netIDs <- dir(x@path, "netID", full.names = TRUE)
  netIDs.success <- file.copy(netIDs, new.ssn.path)
  if(!all(netIDs.success)){
    warning("netID files may not have copied correctly. Try specifying a new ssn path and try again.")
  }
  
  # Then write out the streams dataset as a shapefile
  edges <- as.SpatialLinesDataFrame.SpatialStreamNetwork(x)
  writeOGR(obj = edges, dsn = paste(new.ssn.path, "edges.shp", sep = "/"), layer = "edges", driver = "ESRI Shapefile", ...)
  
  # Output nothing
  invisible()
  
}