#' For each tributary junction in a SpatialStreamNetwork, find the rids of all reaches that meet to form the tributary junction.
#' 
#'@description 
#' 
#'For each tributary junction in a SpatialStreamNetwork, find the rids of all reaches that meet to form the tributary junction.
#' 
#'@usage 
#'    
#'\code{Find.Tribsfunction(ssn.obj, bin.table)} 
#' 
#'@param ssn.obj an object of class SpatialStreamNetwork
#'@param bin.table an object of class data.frame which represents the bindaryID object for a SpatialStreamNetwork
#'@return a data.frame
#'  
#'@export 
Find.Tribs<-function(ssn.obj, bin.table){
  branch1<-numeric()
  branch0<-numeric()
  DSseg<-numeric()
  TribID<-numeric()
  TribUSDist<-numeric()
  tcount<-0
  minShreve<-numeric()
  meanShreve<-numeric()
  
  for(i in 1:nrow(bin.table)){
    b1<-paste(bin.table$binaryID[i], "1", sep="")
    b0<-paste(bin.table$binaryID[i], "0", sep="")
    if(b1 %in% bin.table$binaryID & b0 %in% bin.table$binaryID){
      branch1<-c(branch1, bin.table$rid[bin.table$binaryID==b1])
      branch0<-c(branch0, bin.table$rid[bin.table$binaryID==b0])
      DSseg<-c(DSseg, bin.table$rid[i])
      tcount<-tcount+1
      TribID<-c(TribID, tcount)
      TribUSDist<-c(TribUSDist, ssn.obj@network.line.coords$DistanceUpstream[as.integer(ssn.obj@network.line.coords$SegmentID)==bin.table$rid[i]])
    }
  }
  Tribs<-data.frame(DSseg, branch1,branch0, TribID, TribUSDist)
  Tribs<-Tribs[order(Tribs$TribUSDist),]
  Tribs$TribID<-1:nrow(Tribs)
  return(Tribs)
}