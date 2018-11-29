findOptimalDesignParallel_HPC <- function(
  ssn, 
  new.ssn.path,
  glmssn, 
  afv.column, 
  n.points, 
  utility.function, 
  prior.parameters, 
  n.draws = 500, 
  n.cores = 2,
  extra.arguments = NULL
){
  
  # Check inputs
  
  if(!is(ssn, "SpatialStreamNetwork")){
    stop("The argument 'ssn' must be an object of class SpatialStreamNetwork.")
  }
  n <- nnetwork(ssn)
  if(length(n.points) != 1 & length(n.points) != n){
    stop("n.points must have the same length as the number of networks in ssn, or have a length of 1.")
  }
  if(!is.numeric(n.points) | any(n.points < 1)){
    stop("The argument n.points must have a positive integer value.")
  }
  if(!is.numeric(n.draws) | n.draws != floor(n.draws)){
    stop("The argument n.draws must be a whole-numbered numeric value.")
  }
  if(n.draws < 2){
    stop("Please choose a sensible number of draws for n.draws. Recommended values are 100 or more.")
  }
  ### PARAMETER TO CONTROL WHETHER DESIGNS ARE FOUND BY NETWORK OR ACROSS NETWORKS
  do.separately <- TRUE
  if(length(n.points) == 1){
    do.separately <- FALSE
  }
  if(!is.numeric(n.draws) | n.draws < 1){
    stop("The argument n.draws must have a positive integer value.")
  }
  
  ### Code to detect whether there are replicates on sites
  n.pids <- length(
    unique(
      ssn@obspoints@SSNPoints[[1]]@point.data$pid
    )
  )
  n.locIDs <- length(
    unique(
      ssn@obspoints@SSNPoints[[1]]@point.data$locID
    )
  )
  
  # Extract K for greedy exchange algorithm
  
  if(is.null(extra.arguments)){
    extra.arguments <- list()
  }
  if(is.null(extra.arguments$K)){
    K <- 20 
  } else {
    K <- extra.arguments$K
  }
  
  # Extract important matrices for each network
  
  dist.junc.obs <- getStreamDistMatInOrder(ssn)
  if(length(ssn@predpoints@SSNPoints) > 0){
    dist.junc.prd <- getStreamDistMatInOrder(ssn, "preds")
    dist.junc.pxo <- getStreamDistMatInOrder.predsxobs(ssn)
  }
  
  # Initialise loops by creating a list of potential sites to choose from
  is.replicated <- n.pids != n.locIDs
  if(is.replicated){
    print("Replicates found. Mapping PIDs to locIDs...")
    all.points <- unique(as.character(ssn@obspoints@SSNPoints[[1]]@point.data$locID))
    # Create mapping from locID to pid
    all.locIDs <- as.numeric(as.character(ssn@obspoints@SSNPoints[[1]]@point.data$locID))
    unique.locIDs <- unique(all.locIDs)
    locID.to.pid <- vector("list", n.locIDs)
    for(i in 1:n.locIDs){
      locID.i <- all.locIDs == i
      locID.to.pid[[i]] <- ssn@obspoints@SSNPoints[[1]]@point.data$pid[locID.i]
    }
  } else {
    all.points <- ssn@obspoints@SSNPoints[[1]]@point.data$pid
  }
  
  ### SPLIT FUNCTION HERE INTO TWO CASES:
  ## CASE 1: EACH NETWORK IS TREATED SEPARATELY
  ## CASE 2: ALL NETWORKS ARE TREATED TOGETHER
  
  if(do.separately){
    
    # Create an empty list to store diagnostics
    
    diagnostics <- vector("list", n)
    
    # Prepare to split things by network
    
    network.each.point <- ssn@obspoints@SSNPoints[[1]]@point.data$netID
    if(length(ssn@predpoints@SSNPoints) > 0){
      network.each.pred <- ssn@predpoints@SSNPoints[[1]]@point.data$netID
    }
    
    final.points <- vector(mode = "list", length = n)
    
    # For each network, find the design which optimises the utility function.
    
    for(net in 1:n){
      
      ## Check if the number of points to select in this network is 0
      ## Skip if this is the case
      if(n.points[net] == 0){
        next # Skips entirely
      }
      
      dist.junc.obs.net <- dist.junc.obs[[net]]
      
      if(length(ssn@predpoints@SSNPoints) > 0){
        dist.junc.prd.net <- dist.junc.prd[[net]]
        dist.junc.pxo.a.net <- dist.junc.pxo[[ 2 * net - 1 ]]
        dist.junc.pxo.b.net <- dist.junc.pxo[[ 2 * net ]]
      }
      
      # vector for holding points in this network (locIDs only)
      
      n.final.this.network <- n.points[net]
      final.points.this.network <- vector(mode = "numeric", length = n.final.this.network)
      
      ind1 <- network.each.point == net
      points.this.network <- all.points[ind1]
      points.this.network <- as.numeric(points.this.network)
      n.points.this.network <- length(points.this.network)
      
      ## Check if the network has only one point and only one point is requested.
      ## If true, skip this iteration of the loop
      
      if(n.final.this.network == 1 & n.points.this.network == 1){
        final.points[[net]] <- points.this.network
        next # Keep the single point and skip loop iteration
      }
      
      ## PULL OUT DISTANCE MATRICES HEERE
      ## ONLY DO ONCE PER NETWORK
      
      extra.arguments$Matrices.Obs <- getImportantMatrices.obs(
        dist.junc.obs.net, 
        afv = ssn@obspoints@SSNPoints[[1]]@point.data[ind1,afv.column]
      )
      extra.arguments$net.zero.obs <- matrix(1, n.points.this.network, n.points.this.network)
      
      if(length(ssn@predpoints@SSNPoints) > 0){
        indp <- network.each.pred == net
        extra.arguments$Matrices.prd <- getImportantMatrices.obs(
          dist.junc.prd.net,
          afv = ssn@predpoints@SSNPoints[[1]]@point.data[indp,afv.column]
        )
        np <- length(ssn@predpoints@SSNPoints[[1]]@point.data[indp,afv.column])
        extra.arguments$Matrices.pxo <- getImportantMatrices.pxo(
          dist.junc.pxo.a.net,
          dist.junc.pxo.b.net,
          afv.obs = ssn@obspoints@SSNPoints[[1]]@point.data[ind1,afv.column],
          afv.prd = ssn@predpoints@SSNPoints[[1]]@point.data[indp,afv.column]
        )
        extra.arguments$net.zero.prd <- matrix(1, np, np)
        extra.arguments$net.zero.pxo <- matrix(1, n.points.this.network, np)
      }
      
      # List for diagnostics
      
      diagnostics[[net]] <- vector("list", K)
      
      # vector and list for holding designs in each of K iterations
      
      U.all <- c()
      design.all <- list()
      
      for(k in 1:K){
        
        # Print progress to user
        msg <- paste("Now on hyper-iteration", k, "out of", K)
        print(msg)
        
        # Select random subsample of points
        
        random.sample <- sample(1:n.points.this.network, n.final.this.network, FALSE)
        random.points <- points.this.network[random.sample]
        random.points <- as.character(random.points)
        
        if(is.replicated){
          random.points.to.eval <- c()
          for(i in 1:n.final.this.network){
            locID.ind <- which(unique.locIDs == random.points[i])
            random.points.to.eval <- c(random.points.to.eval, locID.to.pid[[locID.ind]])
            check.ind <- random.points.to.eval %in% row.names(glmssn$sampinfo$X)
            random.points.to.eval <- random.points.to.eval[check.ind]
          }
        } else {
          random.points.to.eval <- random.points
        }
        
        # Evaluate utility once to initialise
        
        U <- utility.function(
          ssn,
          glmssn,
          random.points.to.eval,
          prior.parameters, 
          n.draws,
          extra.arguments
        )
        
        U.ij <- U
        design.now <- list(random.points)
        cond <- max(U.ij) >= U
        last.value <- c()
        designs <- list(random.points)
        
        # Iterate through greedy exchange algorithm
        
        while(cond){
          
          for(i in 1:n.final.this.network){
            
            last.value <- random.points
            ind.i <- !(points.this.network %in% last.value)
            remaining.values <- points.this.network[ind.i]
            remaining.values <- as.character(remaining.values)
            n.j <- length(remaining.values)
            remaining.values <- sample(remaining.values, n.j, FALSE)
            
            ## Precompute the designs to evaluate
            d.list.to.eval <- d.list <- vector("list", n.j)
            for(j in 1:n.j){
              random.points[i] <- remaining.values[j]
              ## Map pids to locIDs if there are temporal replicates
              if(is.replicated){
                random.points.to.eval <- c()
                for(l in 1:n.final.this.network){
                  locID.ind <- which(unique.locIDs == random.points[l])
                  random.points.to.eval <- c(random.points.to.eval, locID.to.pid[[locID.ind]]) # problem
                  check.ind <- random.points.to.eval %in% row.names(glmssn$sampinfo$X)
                  random.points.to.eval <- random.points.to.eval[check.ind]
                  d.list.to.eval[[j]] <- random.points.to.eval
                  d.list[[j]] <- random.points
                }
              } else {
                d.list.to.eval[[j]] <- d.list[[j]] <- random.points
              }
            }
            
            designs <- append(designs, d.list)
            #d.iter <- isplitVector(d.list, chunks = n.cores)
            Uij <- mclapply(
              d.list.to.eval, 
              function(x){
                utility.function(
                  design.points = x,
                  ssn = ssn, 
                  glmssn = glmssn,
                  prior.parameters = prior.parameters,
                  n.draws = n.draws,
                  extra.arguments = extra.arguments
                )
              },
              mc.cores = n.cores
            )
            U.ij <- c(U.ij, unlist(Uij))
            
          }
          # Pump utility values into the list element for this network
          
          design.previous <- design.now
          if(max(U.ij) >= U){
            U <- max(U.ij)
            ind <- U.ij == max(U.ij)
            ind <- (1:length(U.ij))[ind]
            #U.ij <- U.ij[-ind]
            design.now <- tryCatch(
              designs[[ind]], 
              error = function(e){
                warning("Design is not optimal; selected because computation of all utilities failed.")
                return(designs[[2]])
              }
            )
            random.points <- design.now
            
          }
          
          cond <- !all(design.now %in% design.previous)
          
        }
        
        # Print progress to user
        msg <- paste("The best design in this iteration had a utility value of", U)
        print(msg)
        
        ## Store utility values
        
        diagnostics[[net]][[k]] <- U.ij
        
        ## Store max and associated design
        
        U.all <- append(U.all, U)
        design.all <- append(design.all, list(design.now))
        
      }
      
      ind <- (1:length(U.all))[U.all == max(U.all)][1]
      final.points[[net]] <- design.all[[ind]]
    }
    
    ## Get the final set of design points
    final.points <- unlist(final.points) 
    
  } else {
    
    ## One set of diagnostics for each of K hyper-iterations
    
    diagnostics <- vector("list", K)
    
    # Extract out total matrix for the observations
    
    total.mats.obs <- constructTotalMatrix(dist.junc.obs)
    extra.arguments$net.zero.obs <- total.mats.obs$net.zero
    
    # Process matrix to obtain important matrices
    
    extra.arguments$Matrices.Obs <- getImportantMatrices.obs(
      total.mats.obs$d.junc,
      ssn@obspoints@SSNPoints[[1]]@point.data[, afv.column]
    )
    rm(total.mats.obs) # remove from memory in case it is very large
    
    # Do the same for prediction related matrices if needed
    if(length(ssn@predpoints@SSNPoints) > 0){
      total.mats.prd <- constructTotalMatrix(dist.junc.prd)
      extra.arguments$net.zero.prd <- total.mats.prd$net.zero
      extra.arguments$Matrices.prd <- getImportantMatrices.obs(
        total.mats.prd$d.junc,
        ssn@predpoints@SSNPoints[[1]]@point.data[, afv.column]
      )
      rm(total.mats.prd)
      total.mats.pxo <- constructTotalMatrix(dist.junc.pxo, TRUE)
      extra.arguments$net.zero.pxo <- total.mats.pxo$net.zero
      extra.arguments$Matrices.pxo <- getImportantMatrices.pxo(
        total.mats.pxo$d.junc,
        t(total.mats.pxo$d.junc),
        ssn@obspoints@SSNPoints[[1]]@point.data[, afv.column],
        ssn@predpoints@SSNPoints[[1]]@point.data[, afv.column]
      )
      rm(total.mats.pxo)
    }
    
    # Initialise loop
    
    U.all <- c()
    design.all <- list()
    for(k in 1:K){
      
      # Print progress to user
      msg <- paste("Now on hyper-iteration", k, "out of", K)
      print(msg)
      
      # Select random subsample of points
      random.sample <- sample(1:length(all.points), n.points, FALSE)
      random.points <- all.points[random.sample]
      random.points <- as.character(random.points)
      
      if(is.replicated){
        random.points.to.eval <- c()
        for(i in 1:n.points){
          locID.ind <- which(unique.locIDs == random.points[i])
          random.points.to.eval <- c(random.points.to.eval, locID.to.pid[[locID.ind]])
          check.ind <- random.points.to.eval %in% row.names(glmssn$sampinfo$X)
          random.points.to.eval <- random.points.to.eval[check.ind]
        }
      } else {
        random.points.to.eval <- random.points
      }
      
      # Evaluate utility once to initialise
      
      U <- utility.function(
        ssn,
        glmssn,
        random.points.to.eval,
        prior.parameters, 
        n.draws,
        extra.arguments
      )
      
      U.ij <- U
      design.now <- list(random.points)
      cond <- max(U.ij) >= U
      last.value <- c()
      designs <- list(random.points)
      
      while(cond){
        
        for(i in 1:n.points){
          
          last.value <- random.points
          ind.i <- !(all.points %in% last.value)
          remaining.values <- all.points[ind.i]
          remaining.values <- as.character(remaining.values)
          n.j <- length(remaining.values)
          ## Precompute the designs to evaluate
          d.list.to.eval <- d.list <- vector("list", n.j)
          for(j in 1:n.j){
            random.points[i] <- remaining.values[j]
            ## Map pids to locIDs if there are temporal replicates
            if(is.replicated){
              random.points.to.eval <- c()
              for(l in 1:n.points){
                locID.ind <- which(unique.locIDs == random.points[l])
                random.points.to.eval <- c(random.points.to.eval, locID.to.pid[[locID.ind]]) # problem
                check.ind <- random.points.to.eval %in% row.names(glmssn$sampinfo$X)
                random.points.to.eval <- random.points.to.eval[check.ind]
                d.list.to.eval[[j]] <- random.points.to.eval
                d.list[[j]] <- random.points
              }
            } else {
              d.list.to.eval[[j]] <- d.list[[j]] <- random.points
            }
          }
          designs <- append(designs, d.list)
          # d.iter <- isplitVector(d.list, chunks = n.cores)
          Uij <- mclapply(
            d.list.to.eval, 
            function(x){
              utility.function(
                design.points = x,
                ssn = ssn, 
                glmssn = glmssn,
                prior.parameters = prior.parameters,
                n.draws = n.draws,
                extra.arguments = extra.arguments
              )
            },
            mc.cores = n.cores
          )
          U.ij <- c(U.ij, unlist(Uij))
          
        }
        
        design.previous <- design.now
        if(max(U.ij) >= U){
          U <- max(U.ij)
          ind <- U.ij == max(U.ij)
          ind <- (1:length(U.ij))[ind]
          design.now <- tryCatch(
            designs[[ind]], 
            error = function(e){
              warning("Design is not optimal; selected because computation of all utilities failed.")
              return(designs[[2]])
            }
          )
          random.points <- design.now
          
        }
        
        cond <- !all(design.now %in% design.previous)
        
      }
      
      msg <- paste("The best design in this iteration had a utility value of", U)
      print(msg)
      
      # Store diagnostics
      
      diagnostics[[k]] <- U.ij
      
      # Store max and associated design
      
      U.all <- append(U.all, U)
      design.all <- append(design.all, list(design.now))
      
    }
    
    ind.max <- U.all == max(U.all)
    final.points <- unlist(design.all[ind.max])
    
  }
  
  ## create new SSN containing only the selected points
  
  final.points <<- final.points # This needs to be in the global environment for subsetSSN to work? It's really strange.
  # if(is.replicated){
  #   suppressWarnings(subsetSSN(ssn = ssn, filename = new.ssn.path, as.character(locID) %in% final.points))
  # } else {
  #   suppressWarnings(subsetSSN(ssn = ssn, filename = new.ssn.path, pid %in% final.points))
  # }
  # preds <- NULL
  # if(length(ssn@predpoints@SSNPoints) > 0){
  #   preds <- "preds"
  #   ssn.new <- importSSN(new.ssn.path, preds)
  #   createDistMat(ssn.new, preds, TRUE, TRUE)
  # } else{
  #   ssn.new <- importSSN(new.ssn.path, NULL)
  #   createDistMat(ssn.new, NULL, TRUE, FALSE)
  # }
  # rm(final.points, pos = 1) # For some reason, the subsetSSN doesn't work unless final.points is in the global scope. Removing it here.
  
  ### NOTE: ABOVE TEXT IS COMMENTED OUT BECAUSE THE HPC THROWS A FIT WHEN RUNNING SUBSET SSN???
  
  # return updated ssn with other info
  return(
    list(
      # ssn.object.old = ssn,
      # ssn.object.new = ssn.new, 
      networks.separate = do.separately, 
      final.design = final.points,
      utility.values = diagnostics
    )
  )
  
}
