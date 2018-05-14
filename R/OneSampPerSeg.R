#'For a given list of rids, this function randomly selects a single sampling location from each reach.
#' 
#'@description 
#' 
#'For a given list of rids, this function randomly selects a single sampling location from each reach.
#' 
#'@usage 
#'    
#'\code{OneSampPerSeg(ssn.obj, segment.vector)} 
#' 
#'@param ssn.obj an object of class SpatialStreamNetwork
#'@param segment.vector a numeric vector of rids for which samples are desired.
#'@return a numeric vector of pids
#'  
#'@export 
OneSampPerSeg <- function(ssn.obj, segment.vector){
  
  ssn.DF<-ssn.obj@obspoints@SSNPoints[[1]]@point.data
  Sample.Holder <-rep(NA, length(segment.vector))
  
  for (i in 1:length(segment.vector)){
    ssn.DF.now<-subset(ssn.DF, rid==segment.vector[i])
    ## randomly pick a location from segment
    pick1 <- sample(1:nrow(ssn.DF.now), 1)
    Sample.Holder[i] <- ssn.DF.now$pid[pick1]
  } # ends i loop
  return(Sample.Holder)
} 