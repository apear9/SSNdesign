correctSystematicDesign <- function(ssn, spacing, preds = FALSE){
  # Check that number of spacings is either 1 or nnetwork(ssn)
  nn <- nnetwork(ssn)
  nspacing <- length(spacing)
  if(nspacing != 1 & nspacing != nn){
    stop("Either provide one spacing or as many spacings as there are networks in the SSN")
  }
  # Copy spacings if only one spacing given
  if(nspacing == 1){
    spacing <- rep(spacing, nn)
  }
  # get stream distance matrices in order and get hydrological distance matrices
  dm.in.order.obs <- getStreamDistMatInOrder(ssn)
  dh.in.order.obs <- lapply(dm.in.order.obs, function(x){getImportantMatrices.obs(x, rep(1, ncol(x)))$d})
  if(preds){
    dm.in.order.prd <- getStreamDistMatInOrder(ssn, "preds")
    dh.in.order.prd <- lapply(dm.in.order.prd, function(x){getImportantMatrices.obs(x, rep(1, ncol(x)))$d})
  }
  # For each network, find the site with the smallest upstream distance
  dat <- ssn@obspoints@SSNPoints[[1]]@point.data
  smallest.obs <- vector("numeric", nn)
  for(i in 1:nn){
    ind.net <- dat$netID == i
    dat.i <- dat[ind.net, c("pid", "upDist", "netID")]
    smallest.obs[i] <- dat.i$pid[which(dat.i$upDist == min(dat.i$upDist))][1]
  }
  if(preds){
    dat <- ssn@predpoints@SSNPoints[[1]]@point.data
    smallest.prd <- vector("numeric", nn)
    for(i in 1:nn){
      ind.net <- dat$netID == i
      dat.i <- dat[ind.net, c("pid", "upDist", "netID")]
      smallest.prd[i] <- dat.i$pid[which(dat.i$upDist == min(dat.i$upDist))][1]
    }
  }
  # For each network, get sites at specified spacing(s)
  selected.overall.obs <- c()
  for(i in 1:nn){
    best.row <- which(row.names(dh.in.order.obs[[i]]) == smallest.obs[i])
    selected.this.network <- fillSpaceAtIntervals(dh.in.order.obs[[i]], spacing[i], best.row)
    selected.overall.obs <- c(selected.overall.obs, selected.this.network)
  }
  selected.overall.obs <- sort(as.numeric(selected.overall.obs))
  # Same for preds if preds asked for
  if(preds){
    selected.overall.prd <- c()
    for(i in 1:nn){
      best.row <- which(row.names(dh.in.order.prd[[i]]) == smallest.prd[i])
      selected.this.network <- fillSpaceAtIntervals(dh.in.order.prd[[i]], spacing[i], best.row)
      selected.overall.prd <- c(selected.overall.prd, selected.this.network)
    }
    selected.overall.prd <- sort(as.numeric(selected.overall.prd))
  }
  
  # Return in a list
  to.return <- list(
    obs = selected.overall.obs,
    preds = NULL
  )
  if(preds){
    to.return$preds <- selected.overall.prd
  }
  # Return
  return(to.return)
  
}
