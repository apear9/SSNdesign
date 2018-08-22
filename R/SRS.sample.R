#'This function selects a simple random sample of points from a stream network object of class SpatialStreamNetwork.
#' 
#'@description 
#' 
#'This function selects a simple random sample of points from a stream network object of class SpatialStreamNetwork.
#' 
#'@usage 
#'    
#'\code{SRS.sample(ssn.obj, sample.size)} 
#' 
#'@param ssn.obj an object of class SpatialStreamNetwork
#'@param sample.size an integer indicating the desired number of point locations selected for the sample.
#'@return a data.frame
#'  
#'@export 
SRS.sample<-function(ssn.obj, sample.size){
	raw.ssn<-ssn.obj
	raw.df<-raw.ssn@obspoints@SSNPoints[[1]]@point.data
	raw.df$Selected.Sample<-0
	raw.df$Selected.Sample[sort(sample(1:nrow(raw.df),sample.size, replace=FALSE))]<-1
	return(raw.df)
} 