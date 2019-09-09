#' Constructing GRTS designs over time using the 'master sample' approach
#' 
#' @description 
#' 
#' This function takes a \code{SpatialStreamNetwork} object with temporally replicated sites and returns a GRTS design that evolves over time.
#' 
#' The 'master sample' approach was described by Larsen et al. (2008). It exploits the fact that ordered subsets of GRTS designs are themselves GRTS designs in order to sequentially construct GRTS designs.
#' 
#' @param ssn A \code{SpatialStreamNetwork} object
#' @param num.sites A numeric vector. This vector can contain a single element or as many elements as there are temporal steps in the \code{ssn} sites
#' @param rep.variable A string being the variable on which the temporal replicates are defined
#' @return A list which indicates the design for each timestep. The designs are given in lists with two named elements: \code{by.locID}, where the sites are given by their locID values; \code{by.pid}, where the sites are given by their pid values. 
#' 
#' @references 
#' 
#' Larsen, D.P., Olsen, A.R., and Stevens, D.L. (2008). Using a Master Sample to Integrate Stream Monitoring Programs. Journal of Agricultural, Biological and Environmental Statistics, 13(3), 243-254.
#' 
#' @export
evolveGRTSOverTime <- function(ssn, num.sites, rep.variable){
  
  # Input checks
  if(!isSSN(ssn)){
    stop("The argument ssn must be an object of class SpatialStreamNetwork.")
  } # that ssn is of correct class
  ssn.df <- getSSNdata.frame(ssn)
  ssn.locs <- ssn.df$locID
  ssn.pids <- ssn.df$pid
  if(length(unique(ssn.pids)) == length(unique(ssn.locs))){
    stop("There are no temporal replicates of sites in this SSN.")
  } # that ssn has temporal replicates
  
  # Input corrections
  rep.var.lvl <- sort(unique(ssn.df[, rep.variable]))
  if(length(rep.var.lvl) == 1){
    stop("There are no temporal replicates.")
  }
  if(length(num.sites) == 1){
    num.sites <- rep(num.sites[1], length(rep.var.lvl))
  } else if(length(num.sites) == length(rep.var.lvl)) {
    num.sites <- num.sites
  } else {
    stop("The argument num.sites must be length 1 or have the same length as the number of levels in the replication variable.")
  }
  
  # Get total sample size
  samp.sz <- sum(num.sites)
  
  # Begin design process
  dsn.df <- data.frame(
    xcoord = ssn@obspoints@SSNPoints[[1]]@point.coords[,1],
    ycoord = ssn@obspoints@SSNPoints[[1]]@point.coords[,2],
    id = ssn.df$locID,
    locID = ssn.df$locID
  )
  dsn.df <- unique(dsn.df)
  
  design.specs <- list(
    None=list(panel=c(Panel1=samp.sz),seltype="Equal", over=0)
  )
  design.grts <- suppressMessages(suppressWarnings(grts(
    design.specs, 
    DesignID="Site", 
    SiteBegin=1, 
    type.frame="finite",
    src.frame = "att.frame", 
    att.frame = dsn.df, 
    id = "id", 
    xcoord="xcoord",
    ycoord = "ycoord", 
    shapefile=FALSE
  )))
  selected.pool <- anum(design.grts$locID)
  sites.by.timestep <- vector("list", length(num.sites))
  for(i in 1:length(num.sites)){
    sites.by.timestep[[i]] <- list(by.locID = NA, by.pid = NA) 
    num.to.use <- sum(num.sites[1:i])
    if(num.to.use == 0){
      sites.by.timestep[[i]]$by.locID <- numeric(0)
    } else {
      sites.by.timestep[[i]]$by.locID <- selected.pool[1:num.to.use]
    }
  }
  names(sites.by.timestep) <- paste("Period", rep.var.lvl, sep = "_")
  
  # Express designs in terms of pids also
  for(i in 1:length(num.sites)){
    ind.1 <- ssn.df[, rep.variable] == rep.var.lvl[i]
    ind.2 <- ssn.df$locID %in% sites.by.timestep[[i]]$by.locID
    sites.by.timestep[[i]]$by.pid <- ssn.df$pid[ind.1 & ind.2]
    if(i != 1){
      sites.by.timestep[[i]]$by.pid <- c(sites.by.timestep[[i]]$by.pid, sites.by.timestep[[i - 1]]$by.pid)
    }
  }
  
  # Return the design
  return(sites.by.timestep)
  
}