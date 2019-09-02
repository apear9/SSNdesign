#' Find the length of every reach in a SpatialStreamNetwork.
#' 
#'@description 
#' 
#'For a SpatialStreamNetwork, find the length of all reaches.
#'
#'@param raw.ssn an object of class SpatialStreamNetwork
#'@param bin.table an object of class data.frame which represents the bindaryID object for a SpatialStreamNetwork
#'@return a data.frame
#'  
#'@export 
find.segment.lengths<-function(raw.ssn, bin.table){
  
  num.segments<-nrow(raw.ssn@network.line.coords)
  
  ## first: find length of most downstream segment
  st1<-which(raw.ssn@network.line.coords$DistanceUpstream==min(raw.ssn@network.line.coords$DistanceUpstream))
  seg.rid<-raw.ssn@network.line.coords$SegmentID[st1]
  seg.length<-raw.ssn@network.line.coords$DistanceUpstream[st1]
  top.dist<-raw.ssn@network.line.coords$DistanceUpstream[st1]
  bottom.dist<-0
  
  for ( i in (1:nrow(raw.ssn@network.line.coords))[-st1] ){
    seg.rid<-c(seg.rid ,raw.ssn@network.line.coords$SegmentID[i])
    
    rid.below<-find.seg.below(bin.table$binaryID[bin.table$rid==raw.ssn@network.line.coords$SegmentID[i]], bin.table)	
    seg.length<-c(seg.length, (raw.ssn@network.line.coords$DistanceUpstream[i]- 
                                 raw.ssn@network.line.coords$DistanceUpstream[raw.ssn@network.line.coords$SegmentID==rid.below]))
    top.dist<-c(top.dist, raw.ssn@network.line.coords$DistanceUpstream[i])
    bottom.dist<-c(bottom.dist,raw.ssn@network.line.coords$DistanceUpstream[raw.ssn@network.line.coords$SegmentID==rid.below] )
  } 
  rid.length<-data.frame(rid=seg.rid, length=seg.length, TopDist= top.dist, BottomDist=bottom.dist)
  return(rid.length)
} 