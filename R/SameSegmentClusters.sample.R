#'For sampling designs including clusters of samples, this function selects clustered samples within each reach that is provided.
#' 
#'@description 
#' 
#'For sampling designs including clusters of samples, this function selects clustered samples within each reach that is provided.
#' 
#'@usage 
#'    
#'\code{SameSegmentClusters.sample(ssn.obj, ClustDistMethod="prop.shortest.seg", segment.vector, max.dist=NA, cluster.size, start.point.method="random", bin.table)} 
#' 
#'@param ssn.obj an object of class SpatialStreamNetwork
#'@param ClustDistMethod a character vector providing the method for how clustered samples are selected from within a single reach.
#'@param segment.vector a numeric vector of rids for which clustered samples are desired.
#'@param max.dist an optional distance provided as the upper limit for a distance between two locations within a common reach whereby those locations can be selected into a common cluster.
#'@param start.point.method a character vector providing the method for how the first location in each cluster should be selected
#'@param bin.table an object of class data.frame which represents the bindaryID object for a SpatialStreamNetwork
#'@return a data.frame
#'  
#'@export 
SameSegmentClusters.sample<-function(ssn.obj, ClustDistMethod="prop.shortest.seg", segment.vector, max.dist=NA, cluster.size, start.point.method="random", bin.table){
  
  ssn.DF<-ssn.obj@obspoints@SSNPoints[[1]]@point.data
  Sample.Holder<-numeric()
  ssn.DF$Selected.Sample<-0
  
  if(ClustDistMethod=="prop.shortest.seg"){
    # find segment lengths
    seg.lengths<-find.segment.lengths(ssn.obj, bin.table)
    max.dist<-min(seg.lengths$length)
    for (i in 1:length(segment.vector)){
      ssn.DF.now<-subset(ssn.DF, rid==segment.vector[i])
      ## Pick first location
      loc.1<-first.loc(ssn.DF.now, start.point.method)
      
      ## Pick other points in "cluster" of loc.1
      ssn.DF.now$potential.clusterees<-abs(ssn.DF.now$upDist-ssn.DF.now$upDist[ssn.DF.now$pid==loc.1])<=max.dist
      ssn.DF.now$potential.clusterees[ssn.DF.now$pid==loc.1]<-0
      if(cluster.size - 1 > length(ssn.DF.now$pid)) stop("There are not enough sites to include in at least one of the requested clusters. This has most likely happened because your grid of potential sites is too coarse.")
      clust.locs<-sample(ssn.DF.now$pid, (cluster.size-1), prob=ssn.DF.now$potential.clusterees)
      Sample.Holder<-c(Sample.Holder, loc.1, clust.locs)
      
    } # ends i loop
  }
  
  if(ClustDistMethod=="fixed"){
    
    for (i in 1:length(segment.vector)){
      
      ssn.DF.now<-subset(ssn.DF, rid==segment.vector[i])
      ## Pick first location
      loc.1<-first.loc(ssn.DF.now, start.point.method)
      
      ## Pick other points in "cluster" of loc.1
      ssn.DF.now$potential.clusterees<-abs(ssn.DF.now$upDist-ssn.DF.now$upDist[ssn.DF.now$pid==loc.1])<=max.dist
      ssn.DF.now$potential.clusterees[ssn.DF.now$pid==loc.1]<-0
      clust.locs<-sample(ssn.DF.now$pid, (cluster.size-1), prob=ssn.DF.now$potential.clusterees)
      Sample.Holder<-c(Sample.Holder, loc.1, clust.locs)
      
    } # ends i loop
    
  }
  
  if(ClustDistMethod=="NoBound"){
    for (i in 1:length(segment.vector)){
      ssn.DF.now<-subset(ssn.DF, rid==segment.vector[i])
      locs<-sample(ssn.DF.now$pid, cluster.size)
      Sample.Holder<-c(Sample.Holder, locs)
    } # ends i loop			
    
  }
  
  ssn.DF$Selected.Sample[ssn.DF$pid%in%Sample.Holder]<-1
  return(ssn.DF)
} 