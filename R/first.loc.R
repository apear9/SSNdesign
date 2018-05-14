#' Find the first sampling location in a clustered sample
#' 
#'@description 
#' 
#'For a sampling design that includes clusters of sampled locations, this function provides the first sample in the cluster.  
#' 
#'@usage 
#'    
#'\code{first.loc(ssn.DF.now, start.point.method)} 
#' 
#'@param ssn.DF.now an object of class data.frame which is a single reach's subset of the point.data slot fomr a SpatialStreamNetwork 
#'@param start.point.method a character vector indicating which method to apply in selecting the the location, with choices of random, prop.to.end, and prop.to.start
#'@return a numeric scalar which is the selected pid column from the point.data data frame
#'  
#'@export 
first.loc<-function(ssn.DF.now, start.point.method){
  switch(start.point.method,
         random = sample(ssn.DF.now$pid, 1),
         prop.to.end = sample(ssn.DF.now$pid, 1, prob=ssn.DF.now$ratio*((ssn.DF.now$ratio-.5)>0)),
         prop.to.start = sample(ssn.DF.now$pid, 1, prob=(1-ssn.DF.now$ratio)*((ssn.DF.now$ratio-.5)<0))
  ) 
} 