#'Split the temporally replicated point data for a SpatialStreamNetwork object into separate shapefile
#' 
#'@description
#' 
#'Primarily intended to be used with the doAdaptiveDesign function, this function takes an SSN with temporal replicates across some number of sites and splits the sites into separate shapefiles based on the levels of the replication variable. This makes it easier to deal with the SSN in a statistically consistent way when performing sequential design over the network.
#' 
#'@usage 
#' 
#'\code{splitSSNSites(big.ssn, replication.variable, preds = TRUE, preferred.wd = getwd())}
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
#'Good stuff.
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
  
  ## Find number of replicates in the replication variable
  rep.var <- big.ssn@obspoints@SSNPoints[[1]]@point.data[,replication.variable]
  rep.val <- unique(rep.var)
  rep.num <- length(rep.val)
  
  ## Find which value of the replication variable is smallest
  rep.min <- sort(rep.val)[1] # most general, works for factors, character and numeric (though for characters extra care should be taken)
  
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
    }
    prd.ind <- big.dat[, replication.variable] == rep.min
    prd.pid <<- big.dat$pid[prd.ind]
  }
  
  ## Print
  print(paste("Sites have been split by the variable", replication.variable))
  print(paste("New files have been written to the folder", new.wd))
  
  ## Reset wd
  setwd(old.wd)
  
  ## Subset SSN to match the first set of sites
  new.ssn <- subsetSSN(big.ssn, new.ssn.path, subset = pid %in% obs.pid)
  if(preds){
    new.ssn <- subsetPreds(new.ssn, pid %in% prd.pid)
  }
  
  ## Return new ssn object
  
  return(new.ssn)
  
}

