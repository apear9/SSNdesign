constructGRTSOverTime <- function(full.ssn, init.grts, replication.variable, sample.sizes){
  vals <- sort(unique(full.ssn@obspoints@SSNPoints[[1]]@point.data[, replication.variable]))
  reps <- length(vals)
  if(reps == 1){
    stop("No need to use this function to construct a GRTS design. There is only one time step. Use drawStreamNetworkSamples instead.")
  }
  if(length(sample.sizes) == 1){
    sample.sizes <- rep(sample.sizes, reps)
  }
  if(length(sample.sizes) != reps){
    stop("Must specify a single sample size to draw or one for each level of the replication.variable argument.")
  }
  new.full.path <- paste(tempdir(), basename(full.ssn@path), sep = "/")
  writeSSN(full.ssn, new.full.path)
  full.ssn <- updatePath(full.ssn, new.full.path)
  split.path <- paste(tempdir(), "split.ssn", sep = "/")
  preds <- anyPreds(full.ssn)
  split.ssn <- splitSSNSites(full.ssn, split.path, replication.variable, preds, tempdir())
  # Initial GRTS design
  selected <- vector("list", reps)
  if(!missing(init.grts)){
    # Get locIDs already in a GRTS design
    ex.locIDs <<- as.character(
      getSSNdata.frame(init.grts)$locID
    )
    grts.i <- init.grts
  } else {
    new.path <- paste(tempdir(), "grts1.ssn", sep = "/")
    grts.i <- drawStreamNetworkSamples(split.ssn, new.path, TRUE, "GRTS", sample.sizes[1])
    ex.locIDs <<- as.character(
      getSSNdata.frame(grts.i)$locID
    )
  }
  selected[[1]] <- ex.locIDs
  # Rest of the GRTS designs
  for(i in 2:reps){
    # Splice
    new.path <- paste0(tempdir(), "/grts", i, ".ssn")
    grts.i <- spliceSSNSites(grts.i, new.path, paste0("sites", i, ".shp"), paste0("preds", i, ".shp"))
    # Subset to single value of replication variable
    new.path <- paste0(tempdir(), "/grts", i, "select.ssn")
    vals.i <<- vals[i]
    grts.i <- subsetSSN(grts.i, new.path, Year == vals.i & !(locID %in% ex.locIDs))
    # Get grts sample from remaining sites
    new.path <- paste0(tempdir(), "/grts", i, "selected.ssn")
    grts.i <- drawStreamNetworkSamples(grts.i, new.path, TRUE, "GRTS", sample.sizes[i])
    ex.locIDs <<- c(
      ex.locIDs, 
      as.character(
        getSSNdata.frame(grts.i)$locID
      )
    )
    selected[[i]] <- ex.locIDs
  }
  selected <<- selected
  # Subset to fit GRTS designs
  grts.ssns <- vector("list", reps)
  for(i in 1:reps){
    ssn.path <- paste0(goUpOneLevelInSSNPath(full.ssn@path), "grts", i, "design.ssn")
    ind <- rep(FALSE, ndpoints(full.ssn)$inTotal)
    for(j in 1:i){
      ind <- ind | (full.ssn@obspoints@SSNPoints[[1]]@point.data$locID %in% selected[[j]] & full.ssn@obspoints@SSNPoints[[1]]@point.data$Year == vals[j])
    }
    ind <<- ind
    vals.i <<- vals[i]
    grts.ssns[[i]] <- subsetSSN(full.ssn, ssn.path, ind)
    grts.ssns[[i]] <- subsetPreds(grts.ssns[[i]], Year <= vals.i)
  }
  # # Import predpoints if required
  # for(i in 1:reps){
  #   if(i == 1){
  #     if(!missing(init.grts)){
  #       predpts.path <- paste0(init.grts@path, "/preds")
  #     }else{
  #       predpts.path <- paste0(tempdir(), "/grts", 1, ".ssn/preds")
  #     }
  #   } else{
  #     predpts.path <- paste0(tempdir(), "/grts", i, ".ssn/preds")
  #   }
  #   print(predpts.path)
  #   grts.ssns[[i]] <- importPredpts(grts.ssns[[i]], predpts.path, "SSN")
  # }
  # Then return as list
  return(grts.ssns)
}