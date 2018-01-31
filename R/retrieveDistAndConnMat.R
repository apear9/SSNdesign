#' Retrieve distance and connectivity matrices for points in a SpatialStreamNetwork
#' 
#' @\code{retrieveDistAndConnMat} ...
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param networks A numeric or a numeric vector. Defaults to a vector of all the networks in the SpatialStreamNetwork.
#' @return A named list, containing two lists; one contains the distance matrices and the other contains the connectivity matrices.
#' 
#' @examples 
#' \code{#None as yet}
retrieveDistAndConnMat <- function(ssn, networks = sort(as.numeric(as.character(unique((ssn$netID)))))){
  
  # set working directory to "distances/obs" for the ssn in question
  
  oldwd <- getwd()
  tmpwdp <- paste(ssn@path, "distance", "obs", sep = "/")
  if(!dir.exists(tmpwdp)){
    stop("Directory for distance matrices does not exist.")
  }
  tmpwd <- setwd(tmpwdp)
  
  # loop through .Rdata files to create list of distance and connectivity matrices
  n.matrices <- length(networks)
  em <- vector(mode = "list", n.matrices)
  matrices <- list(Distance = em, Connectivity = em)
  for(i in 1:n.matrices){
    netID <- networks[i]
    path.tmp <- paste0("dist.net", netID, ".Rdata")
    if(!file.exists(path.tmp)){
      stop(paste("Distance matrix for network", netID, "does not exist."))
    }
    dm.handle <- file(path.tmp, open = "rb")
    dist.matrix <- unserialize(dm.handle)
    close(dm.handle)
    ordering <- order(as.numeric(rownames(dist.matrix)))
    dist.matrix <- dist.matrix[ordering, ordering]
    matrices$Distance[[i]] <- dist.matrix
    conn.matrix <- 1 - (pmin(dist.matrix, t(dist.matrix)) > 0)* 1
    matrices$Connectivity[[i]] <- conn.matrix
  }
  names(matrices$Distance) <- networks
  names(matrices$Connectivity) <- networks
  # Reset wd and exit
  
  setwd(oldwd)
  return(matrices)
  
}

