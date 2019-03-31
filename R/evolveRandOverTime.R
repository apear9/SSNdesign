#' Constructing random designs over time 
#' 
#' @description 
#' 
#' This function takes a \code{SpatialStreamNetwork} object with temporally replicated sites and returns a random design that evolves over time. Once sites are included in the design, they remain for all subsequent timesteps.
#' 
#' @param ssn A \code{SpatialStreamNetwork} object
#' @param num.sites A numeric vector. This vector can contain a single element or as many elements as there are temporal steps in the \code{ssn} sites
#' @param rep.variable A string being the variable on which the temporal replicates are defined
#' @return A list which indicates the design for each timestep. The designs are given in lists with two named elements: \code{by.locID}, where the sites are given by their locID values; \code{by.pid}, where the sites are given by their pid values. 
#' 
#' @export
evolveRandOverTime <- function(ssn, num.sites, rep.variable){
  
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
  if(samp.sz > length(unique(ssn.df$locID))){
    stop(paste("There aren't enough potential sampling sites to choose from for a sample of size", samp.sz))
  }
  
  # Select sites randomly
  pot.locs <- ssn.df$locID
  pot.locs <- unique(pot.locs)
  # sel.locs <- vector("list", length(rep.var.lvl))
  # for(i in 1:length(rep.var.lvl)){
  #   sel.locs[[i]] <- sample(pot.locs, num.sites[i], FALSE)
  #   pot.locs <- pot.locs[!pot.locs %in% sel.locs[[i]]]
  #   if(i != 1){
  #     sel.locs[[i]] <- c(sel.locs[[i]], sel.locs[[i-1]])
  #   }
  # }
  
  selected.pool <- sample(pot.locs, samp.sz, FALSE)
  sites.by.timestep <- vector("list", length(num.sites))
  for(i in 1:length(num.sites)){
    sites.by.timestep[[i]] <- list(by.locID = NA, by.pid = NA) 
    sites.by.timestep[[i]]$by.locID <- selected.pool[1:sum(num.sites[1:i])]
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