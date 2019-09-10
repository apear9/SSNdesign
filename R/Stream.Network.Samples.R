#'For a SpatialStreamNetwork object with a dense grid of potential sampling locations, this function draws samples according to the sampling designs found most optimal in Som et al. (2014), and a few others.
#' 
#'@description 
#' 
#'For a Spatial Stream Network with a dense grid of potential sampling locations, this function draws samples according to the sampling designs found most optimal in Som et al. (2014), and a few others.
#' 
#'@param ssn.list an object of class list that contains one item names "ssn.object" which is a SpatialStreamNetwork and another item named "bin.table" which is the binaryID table for the assoicated SpatialStreamNetwork
#'@param sample.method a character vector providing the label for the specific sampling design that is desired, with references to Som et al. (2014) including simple random sample "SRS", "GRTS", "GRTSmouth", "GRTSclus", "Headwater.Clust.and.Singles", "Trib.Sets.Head.Singles.sample", "Trib.Sets.Head.Singles.Mouth.sample"
#'@param sample.size a numeric scalar with the desired sample size (USUALLY) but can differ depending on sampling design chosen - something to clean up.
#'@param use.locID a logical indicating whether sites should be selected by locID instead of pid. 
#'@param cluster.number a numeric scalar for the number of sampling locations contained within each clusted sample.
#'@param Wmat leftover from previous version of function when weigths matrix required as input
#'@param DistMat leftover from previous version of function when distance matrix required as input
#'@param ClustDistMethod An argument that controls the way that clusters of sites are created. Defaults to \code{"prop.shortest.seg"}. Other options are \code{"NoBound"} and \code{"fixed"}. If uncertain, leave blank. If you choosed \code{"fixed"}, then you will need to specify \code{max.dist}.
#'@param max.dist The largest distance (in map units) that can separate any two sites belonging to a single cluster. Defaults to the minimum edge length in the SSN. Leave blank if you have no strong beliefs about what this should be. 
#'@param ... Other arguments to internal functions. None implemented as yet.
#'@return a list containing a SpatialStreamNetwork for the selected samples, and a numeric vector with the selected pids
#'  
#'@references 
#' 
#' Som, N.A., Monestiez, P., Ver Hoef, J.M., Zimmerman, D.L., & Peterson, E.E. (2014). Spatial sampling on streams: principles for inference on aquatic networks. \emph{Environmetrics}, \emph{25}(5), 306-323. doi: 10.1002/env.2284.
#'  
#'@export 
Stream.Network.Samples<-function(ssn.list, sample.method, sample.size, use.locID, cluster.number=NULL, Wmat=NULL, DistMat=NULL, ClustDistMethod, max.dist, ...){
  #ssn.object
  ssn.obj<-ssn.list$ssn.object #ssn.object
  bin.table <- do.call(rbind, ssn.list$bin.table) #getBIDtables(ssn.object, networks = 1:nnetwork(ssn.object))
  ssn.obj2<-ssn.obj ## keep an untouched version, potential locs will be dropped from available when moving amonst multiple sample types
  PID.selected<-vector() ## holder of chosen samples through the loops
  ## keep an untouched verson of list too, subset at each sample type
  ssn.list2<-ssn.list
  
  ### keep ssn.list untouched
  ## keep modifying ssn.list2 reducing the selected sample points from future consideration (or selection),
  ## then add full Selected.Sample to ssn.list at end.
  
  for(i in 1:length(sample.method)){
    
    if(sample.method[i]=="SRS"){
      ssn.df2<-SRS.sample(ssn.obj2, sample.size[i])
      PID.selected<-c(PID.selected, ssn.df2$pid[ssn.df2$Selected.Sample==1])
      if(i==1){
        ssn.obj2@obspoints@SSNPoints[[1]]@point.data<-ssn.df2
        tempFile1 <- paste(tempdir(),".ssn",sep="")
        ssn.obj2<-subsetSSN(ssn.obj2, filename=tempFile1,subset=Selected.Sample!= 1)
        ssn.list2$ssn.object<-ssn.obj2
      }
    } # ends SRS
    
    if(sample.method[i]=="GRTS"){
      
      cdat <- ssn.obj2@obspoints@SSNPoints[[1]]@point.data
      cdat$xcoord <- ssn.obj2@obspoints@SSNPoints[[1]]@point.coords[,1] 
      cdat$ycoord <- ssn.obj2@obspoints@SSNPoints[[1]]@point.coords[,2] 
      cdat$id <- cdat$pid
      
      Des1 <- list(None=list(panel=c(Panel1=sample.size[i]),seltype="Equal", over=0))
      Grts2 <- spsurvey::grts(Des1, DesignID="Site", SiteBegin=1, type.frame="finite",
                    src.frame = "att.frame", att.frame = cdat, id = "pid", xcoord="xcoord",
                    ycoord = "ycoord", shapefile=FALSE)
      chosen <- sort(Grts2$id)
      cdat$Selected.Sample <- 0
      cdat$Selected.Sample[cdat$pid %in% chosen] <- 1
      PID.selected<-c(PID.selected, chosen)
      if(i==1){
        ssn.obj2@obspoints@SSNPoints[[1]]@point.data<-cdat
        tempFile1 <- paste(tempdir(),".ssn",sep="")
        ssn.obj2<-subsetSSN(ssn.obj2, filename=tempFile1,subset=Selected.Sample!= 1)
        ssn.list2$ssn.object<-ssn.obj2
      }
    } # ends GRTS
    
    if(sample.method[i]=="GRTSmouth"){
      cdat <- ssn.obj2@obspoints@SSNPoints[[1]]@point.data
      cdat$xcoord <- ssn.obj2@obspoints@SSNPoints[[1]]@point.coords[,1] 
      cdat$ycoord <- ssn.obj2@obspoints@SSNPoints[[1]]@point.coords[,2] 
      cdat$id <- cdat$pid
      
      ## which pid @ mouth?
      which.mouth <- which(cdat$upDist == min(cdat$upDist))
      chosen <- cdat$pid[which.mouth]
      cdat2 <- cdat[-which.mouth,]
      
      sampsize <- sample.size[i] - 1
      Des1 <- list(None=list(panel=c(Panel1=sampsize),seltype="Equal", over=0))
      Grts2 <- spsurvey::grts(Des1, DesignID="Site", SiteBegin=1, type.frame="finite",
                    src.frame = "att.frame", att.frame = cdat2, id = "pid", xcoord="xcoord",
                    ycoord = "ycoord", shapefile=FALSE)
      chosen <- c(chosen, sort(Grts2$id))
      cdat$Selected.Sample <- 0
      cdat$Selected.Sample[cdat$pid %in% chosen] <- 1
      PID.selected<-c(PID.selected, chosen)
      if(i==1){
        ssn.obj2@obspoints@SSNPoints[[1]]@point.data<-cdat
        tempFile1 <- paste(tempdir(),".ssn",sep="")
        ssn.obj2<-subsetSSN(ssn.obj2, filename=tempFile1,subset=Selected.Sample!= 1)
        ssn.list2$ssn.object<-ssn.obj2
      }
    } # ends GRTSmouth
    
    if(sample.method[i]=="GRTSclus"){
      #### cluster.number is the number of number of clusters
      
      cdat <- ssn.obj2@obspoints@SSNPoints[[1]]@point.data
      cdat$xcoord <- ssn.obj2@obspoints@SSNPoints[[1]]@point.coords[,1] #cdat$NEAR_X
      cdat$ycoord <- ssn.obj2@obspoints@SSNPoints[[1]]@point.coords[,2] #cdat$NEAR_Y
      cdat$id <- cdat$pid
      
      if(sample.size[i]%%cluster.number == 0){  ## If sample size divisible by cluster number, then all samples within clusters of the same size
        
        grts.size <- cluster.number
        cluster.size <- sample.size[i]/cluster.number
        Des1 <- list(None=list(panel=c(Panel1=grts.size),seltype="Equal", over=0))
        Grts2 <- spsurvey::grts(Des1, DesignID="Site", SiteBegin=1, type.frame="finite",
                      src.frame = "att.frame", att.frame = cdat, id = "pid", xcoord="xcoord",
                      ycoord = "ycoord", shapefile=FALSE)
        
        ## get the segments from these samples
        Grids <- unique(cdat$rid[cdat$pid %in% Grts2$id])
        ## Get the clusters
        if(missing(ClustDistMethod)) ClustDistMethod <- "prop.shortest.seg"
        if(missing(max.dist)) max.dist <- NA
        cdat<-SameSegmentClusters.sample(ssn.obj2, ClustDistMethod=ClustDistMethod, segment.vector=Grids, 
                                         max.dist=max.dist, cluster.size=cluster.size, start.point.method="random", bin.table=bin.table)
        PID.selected<-c(PID.selected, cdat$pid[cdat$Selected.Sample==1])
        
      }else{  ### here is the case where some GRTS singles and some GRTS clusters
        
        # obtain quantities of samples from GRTS'd clusters or singles
        n.singles <- sample.size[i] - (floor(sample.size[i]/cluster.number)*cluster.number)
        cluster.size <- floor(sample.size[i]/cluster.number)
        grts.size <- cluster.number + n.singles
        Des1 <- list(None=list(panel=c(Panel1=grts.size),seltype="Equal", over=0))
        Grts2 <- spsurvey::grts(Des1, DesignID="Site", SiteBegin=1, type.frame="finite",
                      src.frame = "att.frame", att.frame = cdat, id = "pid", xcoord="xcoord",
                      ycoord = "ycoord", shapefile=FALSE)
        
        ## randomly pull out cluster.number of the chosen pids to cluster with, keep rest as pids
        chosen.pid1 <- sort(Grts2$id)
        forclust <- sample(chosen.pid1, cluster.number, replace=F)
        chosen.pid <- chosen.pid1[!(chosen.pid1 %in% forclust)]
        ## get the segments from these cluster samples
        Grids <- unique(cdat$rid[cdat$pid %in% forclust])
        ## Get the clusters
        if(missing(ClustDistMethod)) ClustDistMethod <- "prop.shortest.seg"
        if(missing(max.dist)) max.dist <- NA
        cdat2<-SameSegmentClusters.sample(ssn.obj2, ClustDistMethod=ClustDistMethod, segment.vector=Grids, 
                                          max.dist=max.dist, cluster.size=cluster.size, start.point.method="random", bin.table=bin.table)
        chosen.pid2 <- cdat2$pid[cdat2$Selected.Sample==1]
        chosen <- c(chosen.pid, chosen.pid2)
        cdat$Selected.Sample <- 0
        cdat$Selected.Sample[cdat$pid %in% chosen] <- 1
        PID.selected<-c(PID.selected, chosen) # Can modify this to achieve the 'by.locID' effect
        
      }
      if(i==1){
        ssn.obj2@obspoints@SSNPoints[[1]]@point.data<-cdat
        tempFile1 <- paste(tempdir(),".ssn",sep="")
        ssn.obj2<-subsetSSN(ssn.obj2, filename=tempFile1,subset=Selected.Sample!= 1)
        ssn.list2$ssn.object<-ssn.obj2
      }
    } # ends GRTSclus
    
    ### Below is option to obtain H1 sample type from the paper
    if(sample.method[i]=="Headwater.Clust.and.Singles"){
      ## will force cluster at mouth segment, other potential clusters from extreme (headwaters) segments, and draw singles from extremes
      ## corrently restricts cluster size to 2 samples per cluster
      
      #Find extreme segments
      ExtSegs<-ssn.list2$ssn.object@data$rid[ssn.list2$ssn.object@data$shreve==1]
      # Mouth = max(shreve)
      mouth.shreve <- max(ssn.list2$ssn.object@data$shreve)
      MouthSeg<-ssn.list2$ssn.object@data$rid[ssn.list2$ssn.object@data$shreve==mouth.shreve]
      
      ## Given sample size, one cluster @ mouth, cluster.number -1 extreme segments, rest of singles from extremes
      extreme.clusters <- sample(ExtSegs, (cluster.number - 1), replace = F)
      
      Clusters <-c(MouthSeg, extreme.clusters)
      
      if(missing(ClustDistMethod)) ClustDistMethod <- "prop.shortest.seg"
      if(missing(max.dist)) max.dist <- NA
      ClustSamps<-SameSegmentClusters.sample(ssn.obj2, ClustDistMethod=ClustDistMethod, segment.vector=Clusters, 
                                             max.dist=max.dist, cluster.size=2, start.point.method="random", bin.table=bin.table)
      chosen.Clust.samps <- ClustSamps$pid[ClustSamps$Selected.Sample ==1]
      
      How.many.samples.left <- sample.size - length(Clusters)*2
      
      ## Need to draw samples, one from all available Extreme segments
      ExtSegsLeft <- ExtSegs[-which(ExtSegs %in% extreme.clusters)]
      Ext.Segs.choose <- sample(ExtSegsLeft, How.many.samples.left, replace=F)
      ExtSamps <- OneSampPerSeg(ssn.obj2, segment.vector = Ext.Segs.choose)
      
      TotalSamps <- c(chosen.Clust.samps, ExtSamps)
      
      PID.selected<-c(PID.selected, TotalSamps)
      if(i==1){
        obj2.DF<-ssn.obj2@obspoints@SSNPoints[[1]]@point.data
        obj2.DF$Selected.Sample <- 0
        obj2.DF$Selected.Sample[obj2.DF$pid %in% TotalSamps] <- 1
        ssn.obj2@obspoints@SSNPoints[[1]]@point.data<-obj2.DF
        tempFile1 <- paste(tempdir(),".ssn",sep="")
        ssn.obj2<-subsetSSN(ssn.obj2, filename=tempFile1,subset=Selected.Sample!= 1)
        ssn.list2$ssn.object<-ssn.obj2
      }
    } # ends if Headwater.Clust.and.Singles
    
    ### Below is option to obtain C3 sample type from the paper
    if(sample.method[i]=="Trib.Sets.Head.Singles.sample"){
      ## Find Tribs
      Trib.DF<-Find.Tribs(ssn.obj2, bin.table)
      ## Select which tribs to sample from
      
      ## For this, start with total.sample size, and cluster.number is those @ tribs, remaining drawn from extreme (headwater) segments, assumes trib clusters contain three samples
      ## cluster number is the number of trib junctions for which cluster samples will be created
      ## pick n/3 tribs
      
      how.many.tribs.2.select<-cluster.number 
      samp.tribs<-sample(1:nrow(Trib.DF), how.many.tribs.2.select)
      tribs.now<-Trib.DF[samp.tribs,]
      ssn.DF<-ssn.obj2@obspoints@SSNPoints[[1]]@point.data
      ssn.DF$Selected.Sample<-0
      PID.here<-numeric()
      RID.here<-numeric()
      for(j in 1:nrow(tribs.now)){
        DS.now<-as.character(tribs.now$DSseg[j])
        DS.dat<-ssn.DF[ssn.DF$rid==DS.now,]
        DS.samp<-DS.dat$pid[DS.dat$ratio==max(DS.dat$ratio)]
        now.1<-as.character(tribs.now$branch1[j])
        dat.1<-ssn.DF[ssn.DF$rid==now.1,]
        samp.1<-dat.1$pid[dat.1$ratio==min(dat.1$ratio)]
        now.0<-as.character(tribs.now$branch0[j])
        dat.0<-ssn.DF[ssn.DF$rid==now.0,]
        samp.0<-dat.0$pid[dat.0$ratio==min(dat.0$ratio)]
        PID.here<-c(PID.here, DS.samp, samp.1, samp.0)
        RID.here<-c(RID.here, as.numeric(as.character(now.1)), as.numeric(as.character(now.0)) )
      }
      
      ## Now, get balance of samples from extremes, not picking extremes that are involved in any trib sets samples
      ## Keep track of branch 1 and branch 0 in a vector.
      
      ExtSegs<-ssn.obj2@data$rid[ssn.obj2@data$shreve==1]
      
      ## Remove trib cluster segs from ExtSegs
      which.in.clust <- which(ExtSegs %in% RID.here)
      ExtSegs2 <- ExtSegs[-which.in.clust]
      
      samples.left <- sample.size - 3*cluster.number
      
      ## Randomly pick samples.left number of extreme segments
      Ext.SingleSamp.Segs <- sample(ExtSegs2, samples.left)
      
      Ext.SingleSamps <- OneSampPerSeg(ssn.obj2, segment.vector = Ext.SingleSamp.Segs)
      
      TotalSamps <- c(PID.here, Ext.SingleSamps)
      
      PID.selected<-c(PID.selected, TotalSamps)
      
      if(i==1){
        obj2.DF<-ssn.obj2@obspoints@SSNPoints[[1]]@point.data
        obj2.DF$Selected.Sample <- 0
        obj2.DF$Selected.Sample[obj2.DF$pid %in% TotalSamps] <- 1
        ssn.obj2@obspoints@SSNPoints[[1]]@point.data<-obj2.DF
        tempFile1 <- paste(tempdir(),".ssn",sep="")
        ssn.obj2<-subsetSSN(ssn.obj2, filename=tempFile1,subset=Selected.Sample!= 1)
        ssn.list2$ssn.object<-ssn.obj2
      }
      
    }# ends Trib.Sets.Head.Singles
    
    
    ### Below is option to obtain C4 sample type from the paper    
    if(sample.method[i]=="Trib.Sets.Head.Singles.Mouth.sample"){
      ## Find Tribs
      Trib.DF<-Find.Tribs(ssn.obj2, bin.table)
      
      ## Select which tribs to sample from
      ## Start with total.sample size, and cluster.number is those at tribs, with what's leftover in extreme (headwater) segments, make sure divisible by 3 here
      how.many.tribs.2.select<-cluster.number 
      samp.tribs<-sample(1:nrow(Trib.DF), how.many.tribs.2.select)
      tribs.now<-Trib.DF[samp.tribs,]
      #ssn.DF<-getSSNdata.frame(ssn.obj2)
      ssn.DF<-ssn.obj2@obspoints@SSNPoints[[1]]@point.data
      ssn.DF$Selected.Sample<-0
      PID.here<-numeric()
      RID.here<-numeric()
      for(j in 1:nrow(tribs.now)){
        DS.now<-as.character(tribs.now$DSseg[j])
        DS.dat<-ssn.DF[ssn.DF$rid==DS.now,]
        DS.samp<-DS.dat$pid[DS.dat$ratio==max(DS.dat$ratio)]
        now.1<-as.character(tribs.now$branch1[j])
        dat.1<-ssn.DF[ssn.DF$rid==now.1,]
        samp.1<-dat.1$pid[dat.1$ratio==min(dat.1$ratio)]
        now.0<-as.character(tribs.now$branch0[j])
        dat.0<-ssn.DF[ssn.DF$rid==now.0,]
        samp.0<-dat.0$pid[dat.0$ratio==min(dat.0$ratio)]
        PID.here<-c(PID.here, DS.samp, samp.1, samp.0)
        RID.here<-c(RID.here, as.numeric(as.character(now.1)), as.numeric(as.character(now.0)) )
      }
      
      ## Now, get balance of samples from extremes, not picking extremes that are involved in any trib sets samples
      ## Keep track of branch 1 and branch 0 in a vector.
      
      ExtSegs<-ssn.obj2@data$rid[ssn.obj2@data$shreve==1]
      
      ## Remove trib cluster segs from ExtSegs
      which.in.clust <- which(ExtSegs %in% RID.here)
      ExtSegs2 <- ExtSegs[-which.in.clust]
      
      samples.left <- sample.size - 3*cluster.number - 1
      
      ## Randomly pick samples.left number of extreme segments
      Ext.SingleSamp.Segs <- sample(ExtSegs2, samples.left)
      
      Ext.SingleSamps <- OneSampPerSeg(ssn.obj2, segment.vector = Ext.SingleSamp.Segs)
      
      TotalSamps <- c(PID.here, Ext.SingleSamps)
      
      ## Add location at mouth (lowest sample point in network among all points)
      which.mouth <- which(ssn.obj2@obspoints@SSNPoints[[1]]@point.data$upDist== min(ssn.obj2@obspoints@SSNPoints[[1]]@point.data$upDist))
      mouth.pid <- ssn.obj2@obspoints@SSNPoints[[1]]@point.data$pid[which.mouth]
      
      TotalSamps <-c(mouth.pid, TotalSamps)
      
      PID.selected<-c(PID.selected, TotalSamps)
      
      if(i==1){
        obj2.DF<-ssn.obj2@obspoints@SSNPoints[[1]]@point.data
        obj2.DF$Selected.Sample <- 0
        obj2.DF$Selected.Sample[obj2.DF$pid %in% TotalSamps] <- 1
        ssn.obj2@obspoints@SSNPoints[[1]]@point.data<-obj2.DF
        tempFile1 <- paste(tempdir(),".ssn",sep="")
        ssn.obj2<-subsetSSN(ssn.obj2, filename=tempFile1,subset=Selected.Sample!= 1)
        ssn.list2$ssn.object<-ssn.obj2
      }
      
    }# ends Trib.Sets.Head.SinglesMouth
    
  } # ends i loop
  
  ## Here is where to put all samples selected from each [i] back together
  
  ssn.obj@obspoints@SSNPoints[[1]]@point.data$Selected.Sample<-0
  ssn.obj@obspoints@SSNPoints[[1]]@point.data$Selected.Sample[ssn.obj@obspoints@SSNPoints[[1]]@point.data$pid %in% PID.selected]<-1
  tempFile1 <- paste(tempdir(),".ssn",sep="")
  ssn.obj<-subsetSSN(ssn.obj, filename=tempFile1, subset = Selected.Sample == 1)
  list(ssn.obj.samples = ssn.obj, Selected.PID = PID.selected)
} 
