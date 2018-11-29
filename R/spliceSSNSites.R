#'Introduce new sites into a SpatialStreamNetwork object from a shapefile
#' 
#'@description 
#' 
#' This function does the reverse of \code{\link{splitSSNSites}}. Instead of separating a \code{SpatialStreamNetwork} object into different shapefiles, this function combines a set of observed points from a shapefile with an existing SSN object.
#' 
#'@usage 
#' 
#' \code{spliceSSNSites(ssn, splice.obs, splice.preds = NULL)}
#' 
#'@param ssn An object of class SpatialStreamNetwork
#'@param new.ssn.path A file path to a folder to store the augmented ssn.
#'@param splice.obs A file path for the shapefile of observed sites which should be brought into the ssn object.
#'@param splice.preds Optionally a file path for a shapefile of prediction points which should also be brought into the ssn object. This can be (and is by default) NULL if there are no such prediction points.
#'@return An object of class SpatialStreamNetwork that contains the old and new sites.
#' 
#'@details 
#' 
#' As for \code{splitSSNSites}, this function is wrapped by \code{\link{doAdaptiveDesign}}. While \code{doAdaptiveDesign} should be sufficient for solving most adaptive design problems, more specialised or unusual applicaitons may require this function to be used separately.  
#'  
#'@examples
#'
#'\dontrun{#code}
#'  
#'@export
spliceSSNSites <- function(
  ssn, 
  new.ssn.path,
  splice.obs,
  splice.preds = NULL 
){
  
  ## Check that inputs are correct
  
  # ssn must be a spatialstreamnetwork
  if(class(ssn) != "SpatialStreamNetwork"){
    stop("The argument ssn must be a SpatialStreamNetwork")
  }
  # Files must exist if specified
  if(!file.exists(splice.obs)){
    stop("Observed sites shapefile does not exist!")
  }
  preds <- FALSE
  if(!is.null(splice.preds)){
    preds <- TRUE
    if(!file.exists(splice.preds)){
      stop("Prediction points shapefile does not exist!")
    }
  }
  
  ## Create new directory
  if(!dir.exists(new.ssn.path)){
    dir.create(new.ssn.path)
  }
  
  ## Get proj4string from ssn
  p4s.old <- ssn@proj4string
  
  ## Ingest the shapefile of new observed (plus maybe prediction) sites 
  lyr.obs <- tools::file_path_sans_ext(basename(splice.obs))
  shp.obs <- readOGR(dsn = splice.obs, layer = lyr.obs)
  shp.obs <- spTransform(shp.obs, p4s.old) # project in case this has different CRS
  # overwrite row names with pids
  row.names(shp.obs@data) <- row.names(shp.obs@coords) <- as.character(shp.obs@data$pid)
  if(preds){
    lyr.prd <- tools::file_path_sans_ext(basename(splice.preds))
    shp.prd <- readOGR(dsn = splice.preds, layer = lyr.prd)
    shp.prd <- spTransform(shp.prd, p4s.old) # as above
    # overwrite row names with pids
    row.names(shp.prd@data) <- row.names(shp.prd@coords) <- as.character(shp.prd@data$pid)
  }
  
  ## Ingest sites of existing ssn as shapefile
  dsn.obs.old <- paste0(ssn@path, "/sites.shp")
  shp.obs.old <- readOGR(dsn.obs.old, "sites")
  # overwrite row names with pids
  row.names(shp.obs.old@data) <- row.names(shp.obs.old@coords) <- as.character(shp.obs.old@data$pid)
  if(preds){
    dsn.prd.old <- paste0(ssn@path, "/preds.shp")
    shp.prd.old <- readOGR(dsn.prd.old, "preds")
    # overwrite row names with pids
    row.names(shp.prd.old@data) <- row.names(shp.prd.old@coords) <- as.character(shp.prd.old@data$pid)
  }
  
  ## Recombine existing sites plus new ones (inside SSN)
  dat.old <- shp.obs.old@data
  crd.old <- shp.obs.old@coords
  nam.old <- names(dat.old)
  dat.new <- shp.obs@data
  crd.new <- shp.obs@coords
  nam.new <- names(shp.obs)
  
  # Check that all names are shared
  frwd <- all(nam.new %in% nam.old)
  bkwd <- all(nam.old %in% nam.new)
  if(!frwd){
    stop("Column names of new and old data for the prediction sites do not match. Please check this and adjust the names in the shapefile accordingly.")
  }
  if(!bkwd){
    dat.old <- dat.old[, nam.old %in% nam.new]
  }
  
  ## Create new sites SpatialPointsDataFrame
  new.dat.obs <- rbind(dat.old, dat.new)
  new.crd.obs <- rbind(crd.old, crd.new)
  pid.odr <- match(sort(new.dat.obs$pid), new.dat.obs$pid)
  new.dat.obs <- new.dat.obs[pid.odr, ]
  new.crd.obs <- new.crd.obs[pid.odr, ]
  obs.new <- SpatialPointsDataFrame(new.crd.obs, new.dat.obs, proj4string = p4s.old)
  
  ## Repeat process for preds if necessary
  if(preds){
    ## Recombine existing sites plus new ones (inside SSN)
    dat.old <- shp.prd.old@data
    crd.old <- shp.prd.old@coords
    nam.old <- names(dat.old)
    dat.new <- shp.prd@data
    crd.new <- shp.prd@coords
    nam.new <- names(shp.prd)
    
    # Check that all names are shared
    frwd <- all(nam.new %in% nam.old)
    bkwd <- all(nam.old %in% nam.new)
    if(!frwd){
      stop("Column names of new and old data for the prediction sites do not match. Please check this and adjust the names in the shapefile accordingly.")
    }
    if(!bkwd){
      dat.old <- dat.old[, nam.old %in% nam.new]
    }
    
    ## Create new sites SpatialPointsDataFrame
    new.dat.prd <- rbind(dat.old, dat.new)
    new.crd.prd <- rbind(crd.old, crd.new)
    pid.odr <- match(sort(new.dat.prd$pid), new.dat.prd$pid)
    new.dat.prd <- new.dat.prd[pid.odr, ]
    new.crd.prd <- new.crd.prd[pid.odr, ]
    prd.new <- SpatialPointsDataFrame(new.crd.prd, new.dat.prd, proj4string = p4s.old)
    
  }
  
  ## Copy old files
  files.to.copy <- dir(ssn@path, full.names = TRUE) # list full file paths
  file.copy(files.to.copy, new.ssn.path) # by default, folders not copied so distance matrices will remain in old directory
  
  ## Shoehorn the new sites shapefile over the top of the old one
  shp.name <- paste0(new.ssn.path, "/sites.shp")
  writeOGR(obs.new, shp.name, "sites", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  SSN:::write.dbf.SSN(new.dat.obs, "sites", max_nchar = 30)
  
  ## Same for predpoints if they exist
  if(preds){
    shp.name <- paste0(new.ssn.path, "/preds.shp")
    writeOGR(prd.new, shp.name, "preds", driver = "ESRI Shapefile", overwrite_layer = TRUE)
    SSN:::write.dbf.SSN(new.dat.prd, "preds", max_nchar = 30)
  }
  
  ## Import the SSN for hard reset of pid ordering, etc.
  if(preds){
    predpts <- "preds"
  } else {
    predpts <- NULL
  }
  new.ssn <- importSSN(new.ssn.path,predpts)
  
  ## Spit out updated ssn
  return(new.ssn)

}
