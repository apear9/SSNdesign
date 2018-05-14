#' Find the rid for the reach immediatley downstream of another reach.
#' 
#'@description 
#' 
#'For a reach in a SpatialStreamNetwork, this function finds the rid for the reach immediately downstream.
#' 
#'@usage 
#'    
#'\code{find.seg.below(x, IDtable)} 
#' 
#'@param x a numeric scaler which represents the binaryID address for a rich object of class data.frame which is a single reach's subset of the point.data slot fomr a SpatialStreamNetwork 
#'@param IDtable an object of class data.frame which represents the bindaryID object for a SpatialStreamNetwork
#'@return a numeric scalar which is the rid for the reach immediatley downstream of the reach whose address is provided as x
#'  
#'@export 
find.seg.below<-function(x, IDtable){
  # take off last character of x to find seg directly below
  seg.below<-substr(x, 1, nchar(x)-1)
  # get rid of seg below
  rid.below<-IDtable$rid[which(IDtable$binaryID==seg.below)]
  return(rid.below)
} 