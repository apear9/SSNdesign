#'Split the temporally replicated point data for a SpatialStreamNetwork object into separate shapefile
#' 
#'@description
#' 
#'Primarily intended to be used with the \code{\link{optimiseSSNDesign}} function, this function takes an SSN with temporal replicates across some number of sites and splits the sites into separate shapefiles based on the levels of the replication variable. This makes it easier to deal with the SSN in a statistically consistent way when performing sequential design over the network.
#' 
#'@param big.ssn The SpatialStreamNetwork object which contains the temporally replicated site data.
#'@param new.ssn.path A file path for a new SSN object.
#'@param replication.variable A string which indicates which variable in the \code{big.ssn@obspoints@SSNPoints[[1]]@point.data} should be used to split big.ssn into separate shapefiles.
#'@param preds A logical which indicates whether the prediction sites should be dealt with in the same way. This predpoints object should contain the same levels of the replication variable as the obspoints slot.
#'@param preferred.wd A directory where the shapefiles should be written.
#'@return This function returns an object of class SpatialStreamNetwork. This SpatialStreamNetwork will contain the sites with the lowest value of the replication variable.
#' 
#'@details 
#' 
#'This function is wrapped by \code{\link{optimiseSSNDesign}}. The \code{optimiseSSNDesign} function should be sufficient for most adaptive design problems where this is required. However, for more bespoke applications, it may be worthwhile to use this function on its own instead. 
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#' s <- createSSN(10, binomialDesign(100, 2), path = tempPath("r.ssn"), importToR = T)
#' 
#' # Split SSN sites into two shapefiles, one for each year
#' split <- splitSSNSites(s, tempPath("split.ssn"), "Time", FALSE, tempdir())
#' 
#' # Join the year 2 shapefile back to the year 1 shapefile in the SSN
#' spliced <- spliceSSNSites(split, tempPath("spliced.ssn"), paste0(tempdir(), "/sites2.shp"))
#' 
#' }
#' 
#'@export
splitSSNSites <- function(
  big.ssn,
  new.ssn.path,
  replication.variable,
  preds = TRUE,
  preferred.wd = getwd()
){
  
  ## Check that inputs are correct
  
  # ssn must be a spatialstreamnetwork
  if(class(big.ssn) != "SpatialStreamNetwork"){
    stop("The argument ssn must be a SpatialStreamNetwork")
  }
  
  if(!dir.exists(preferred.wd)){
    stop("Please specify a valid directory.")
  }
  
  ## Save old working directory before setting it to new working directory
  old.wd <- getwd()
  new.wd <- preferred.wd
  setwd(new.wd)
  ## Reset wd on exit
  on.exit(
    setwd(old.wd)
  )
  
  ## Find number of replicates in the replication variable
  rep.var <- big.ssn@obspoints@SSNPoints[[1]]@point.data[,replication.variable]
  rep.val <- unique(rep.var)
  rep.num <- length(rep.val)
  
  ## Find which value of the replication variable is smallest
  rep.min <- sort(rep.val)[1] # most general, works for factors, character and numeric
  
  ## Start splitting (observations only at this stage)
  prj.4st <- big.ssn@proj4string
  big.dat <- big.ssn@obspoints@SSNPoints[[1]]@point.data
  big.crd <- big.ssn@obspoints@SSNPoints[[1]]@point.coords
  for(i in 1:rep.num){
    ind.i <- big.dat[, replication.variable] == rep.val[i]
    dat.i <- big.dat[ind.i, ]
    pid.i <- dat.i$pid
    crd.i <- big.crd[row.names(big.crd) %in% pid.i,] # DO NOT TRUST THE ORDERING OF THE COORDINATES THESE THINGS ARE EVIL
    spd.i <- SpatialPointsDataFrame(coords = crd.i, data = dat.i, proj4string = prj.4st)
    writeOGR(spd.i, paste0("sites", i, ".shp"), paste0("sites", i), "ESRI Shapefile", overwrite_layer = TRUE)
  }
  obs.ind <- big.dat[, replication.variable] == rep.min
  obs.pid <<- big.dat$pid[obs.ind]
  
  ## Subset SSN to match the first set of sites
  new.ssn <- suppressWarnings(suppressMessages(subsetSSN(big.ssn, new.ssn.path, subset = pid %in% obs.pid)))
  rm(
    list = "obs.pid", 
    pos = ".GlobalEnv"
  )
  
  ## Do same for preds if they exist
  if(preds){
    prj.4st <- big.ssn@proj4string
    big.dat <- big.ssn@predpoints@SSNPoints[[1]]@point.data
    big.crd <- big.ssn@predpoints@SSNPoints[[1]]@point.coords
    for(i in 1:rep.num){
      ind.i <- big.dat[, replication.variable] == rep.val[i]
      dat.i <- big.dat[ind.i, ]
      pid.i <- dat.i$pid
      crd.i <- big.crd[row.names(big.crd) %in% pid.i,] # DO NOT TRUST THE ORDERING OF THE COORDINATES THESE THINGS ARE EVIL
      spd.i <- SpatialPointsDataFrame(coords = crd.i, data = dat.i, proj4string = prj.4st)
      writeOGR(spd.i, paste0("preds", i, ".shp"), paste0("preds", i), "ESRI Shapefile", overwrite_layer = TRUE)
      if(i == 1){
        writeOGR(spd.i, paste(new.ssn@path, "preds.shp", sep = "/"), "preds", "ESRI Shapefile", overwrite_layer = TRUE)
        new.ssn <- suppressWarnings(suppressMessages(importSSN(new.ssn@path, "preds")))
      }
    }
  }
  
  ## Print messages to user
  message(paste("Sites have been split by the variable", replication.variable))
  message(paste("New files have been written to the folder", new.wd))
  
  ## Return new ssn object
  return(new.ssn)
  
}

