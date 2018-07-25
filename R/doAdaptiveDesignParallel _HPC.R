#' A function to solve adaptive / sequential design problems on SpatialStreamNetworks
#'@description 
#'
#'Finds an optimal design consisting of a specified number of points on a SpatialStreamNetwork, given a set of sites that must remain fixed. 
#' 
#'@usage 
#' 
#'\code{doAdaptiveDesignParallel(ssn, glmssn, fixed.points, afv.column, n.points, utility.function, prior.parameters, n.draws = 500, n.cores = 2, extra.arguments = NULL)}
#' 
#'@param ssn an object of class SpatialStreamNetwork
#'@param new.ssn.path a character string specifying the ssn folder in which to store the output
#'@param glmssn an object of class glmssn
#'@param fixed.points a numeric vector of the pids for sites which must remain fixed in the design
#'@param afv.column the name of the column in the SpatialStreamNetwork object that contains the additive function values
#'@param n.points the number of points to be included in the final design. This can be a single number, in which case all networks will have the same number of points. This can also be a vector with the same length as the number of networks in ssn. 
#'@param utility.function a function with the signature utility.function. Users may define their own. This package provides several built-in utility functions, including \code{\link{sequentialDOptimality}}, \code{\link{KOptimality}}, and \code{\link{sequentialCPOptimality}}. 
#'@param prior.parameters a function to act as a prior for covariance parameter values
#'@param n.draws a numeric value for the number of Monte Carlo draws to take when evaluating potential designs
#'@param n.cores a numeric value indicating the number of cores on which to parallelise the process.
#'@param extra.arguments a list of any extra parameters which can be used to control the behaviour of this function or the utility function
#'@return An object of class SpatialStreamNetwork. The SSNPoints for the obspoints slot will be updated to reflect the selected design. 
#'
#'@details 
#'
#'This function implements the greedy exchange algorithm per Falk et al. (2014) to quickly find an optimal design from a set of sites. The algorithm does not always discover the absolute best design; however, even when this occurs, the algorithm will quickly discover a reasonably effective design. 
#'  
#'@export
doAdaptiveDesignParallel_HPC <- function(
  ssn,
  new.ssn.path,
  glmssn,
  fixed.points,
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
  
  # Overwrite model matrix in the model with the fuller model matrix for the SSN (necessary for utility functions)
  
  datXY <- SSN:::dataXY(
    glmssn$args$formula,
    ssn@obspoints@SSNPoints[[1]]@point.data,
    glmssn$args$family,
    glmssn$args$trialscol,
    glmssn$args$trans.power,
    glmssn$args$trans.shift,
    order(ssn@obspoints@SSNPoints[[1]]@point.data$pid),
    glmssn$args$CorModels
  )
  
  glmssn$sampinfo$X <- datXY$Xmats$X1
  
  # Extract important matrices for each network
  
  dist.junc.obs <- getStreamDistMatInOrder(ssn)
  if(length(ssn@predpoints@SSNPoints) > 0){
    dist.junc.prd <- getStreamDistMatInOrder(ssn, "preds")
    dist.junc.pxo <- getStreamDistMatInOrder.predsxobs(ssn)
  }
  
  # Initialise for loop
  
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
      locID.i <- all.locIDs == unique.locIDs[i]
      locID.to.pid[[i]] <- ssn@obspoints@SSNPoints[[1]]@point.data$pid[locID.i]
    }
    all.points <- all.points[!(all.points %in% fixed.points)]
  } else {
    all.points <- ssn@obspoints@SSNPoints[[1]]@point.data$pid
    all.points <- all.points[!(all.points %in% fixed.points)]
  }
  
  #### SPLIT HERE TO EITHER RUN EACH NETWORK SEPARATELY OR RUN ALL SITES TOGETHER
  
  if(do.separately){
    
    # Create an empty list to store diagnostics
    
    diagnostics <- vector("list", n)
    
    # Prepare to split stuff by network
    
    network.each.point <- ssn@obspoints@SSNPoints[[1]]@point.data$netID
    if(length(ssn@predpoints@SSNPoints) > 0){
      network.each.pred <- ssn@predpoints@SSNPoints[[1]]@point.data$netID
    }
    final.points <- vector(mode = "list", length = n)
    
    # For each network, find the design which optimises the utility function.
    for(net in 1:n){
      
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
      ind2 <- !(points.this.network %in% fixed.points)
      fixed.points.this.network <- points.this.network[!ind2]
      extra.arguments$fixed <- fixed.points.this.network
      points.this.network <- points.this.network[ind2]
      n.points.this.network <- length(points.this.network)
      
      # Create net.zero.obs matrix, with as many rows and columns as sites in the SSN object for the ith network
      nrows <- sum(ind1)
      extra.arguments$net.zero.obs <- matrix(1, nrows, nrows)
      
      ## PULL OUT DISTANCE MATRICES HEERE
      ## ONLY DO ONCE PER NETWORK
      
      extra.arguments$Matrices.Obs <- getImportantMatrices.obs(
        dist.junc.obs.net, 
        afv = ssn@obspoints@SSNPoints[[1]]@point.data[ind1,afv.column]
      )
      
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
      
      if(is.replicated){
        # Map fixed points locIDs to pids
        fixed.points.to.eval <- c()
        for(i in 1:length(fixed.points.this.network)){
          locID.ind <- which(unique.locIDs == fixed.points.this.network[i])
          fixed.points.to.eval <- c(fixed.points.to.eval, locID.to.pid[[locID.ind]])
          extra.arguments$fixed <- fixed.points.to.eval
        }
      } else {
        extra.arguments$fixed <- fixed.points.this.network
      }
      
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
        
        # Join locIDs back to PIDs if req'd
        if(is.replicated){
          random.points.to.eval <- c()
          for(i in 1:n.final.this.network){
            locID.ind <- which(unique.locIDs == random.points[i])
            random.points.to.eval <- c(random.points.to.eval, locID.to.pid[[locID.ind]])
          }
          check.ind <- random.points.to.eval %in% row.names(glmssn$sampinfo$X)
          random.points.to.eval <- random.points.to.eval[check.ind]
        } else {
          random.points.to.eval <- random.points
        }
        random.and.fixed.points <- c(random.points, as.character(fixed.points.this.network))
        random.and.fixed.points.to.eval <- c(random.points.to.eval, extra.arguments$fixed)
        
        # Evaluate utility once to initialise
        U <- utility.function(
          ssn,
          glmssn,
          random.and.fixed.points.to.eval,
          prior.parameters, 
          n.draws,
          extra.arguments
        )
        
        U.ij <- U
        design.now <- list(random.and.fixed.points)
        cond <- max(U.ij) >= U
        last.value <- c()
        designs <- list(random.and.fixed.points)
        
        # Iterate through greedy exchange algorithm
        
        while(cond){
          
          for(i in 1:n.final.this.network){
            
            last.value <- random.points
            ind.i <- !(points.this.network %in% last.value)
            remaining.values <- points.this.network[ind.i]
            remaining.values <- as.character(remaining.values)
            n.j <- length(remaining.values)

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
                  d.list.to.eval[[j]] <- c(random.points.to.eval, extra.arguments$fixed)
                  d.list[[j]] <- c(random.points, fixed.points.this.network)
                }
              } else {
                d.list.to.eval[[j]] <- d.list[[j]] <- c(random.points, extra.arguments$fixed)
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
          #cnt = cnt + 1
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
    
    ## create new SSN containing only the selected points
    
    final.points <- unlist(final.points)
    
  } else {
    
    ## Define one set of diagnostics for each of K hyper-iterations
    
    diagnostics <- vector("list", K)
    
    # Set fixed points inside extra.arguments for utility functions to retrieve if necessary
    
    if(is.replicated){
      # Map fixed points locIDs to pids
      fixed.points.to.eval <- c()
      for(i in 1:length(fixed.points)){
        locID.ind <- which(unique.locIDs == fixed.points[i])
        fixed.points.to.eval <- c(fixed.points.to.eval, locID.to.pid[[locID.ind]])
        extra.arguments$fixed <- fixed.points.to.eval
      }
    } else {
      extra.arguments$fixed <- fixed.points
    }
    
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
      random.points <- sample(all.points, n.points, FALSE)
      #random.points <- all.points[random.sample]
      # Join locIDs back to PIDs if req'd
      if(is.replicated){
        random.points.to.eval <- c()
        for(i in 1:n.points){
          locID.ind <- which(unique.locIDs == random.points[i])
          random.points.to.eval <- c(random.points.to.eval, locID.to.pid[[locID.ind]])
        }
        check.ind <- random.points.to.eval %in% row.names(glmssn$sampinfo$X)
        random.points.to.eval <- random.points.to.eval[check.ind]
      } else {
        random.points.to.eval <- random.points
      }
      random.and.fixed.points <- c(random.points, as.character(fixed.points))
      random.and.fixed.points.to.eval <- c(random.points.to.eval, extra.arguments$fixed)
      
      # Evaluate utility once to initialise
      
      U <- utility.function(
        ssn,
        glmssn,
        random.and.fixed.points.to.eval,
        prior.parameters, 
        n.draws,
        extra.arguments
      )
      
      U.ij <- U
      design.now <- list(random.and.fixed.points)
      cond <- max(U.ij) >= U
      last.value <- c()
      designs <- list(random.and.fixed.points)
      
      # Iterate through greedy exchange algorithm
      
      while(cond){
        
        for(i in 1:n.points){
          
          last.value <- random.points
          ind.i <- !(all.points %in% last.value)
          remaining.values <- all.points[ind.i]
          remaining.values <- as.character(remaining.values)
          n.j <- length(remaining.values)
          
          d.list <- d.list.to.eval <- vector("list", n.j)
          for(j in 1:n.j){
            random.points[i] <- remaining.values[j]
            ## Map pids to locIDs if there are temporal replicates
            if(is.replicated){
              random.points.to.eval <- c()
              for(l in 1:n.points){
                locID.ind <- which(unique.locIDs == random.points[l])
                random.points.to.eval <- c(random.points.to.eval, locID.to.pid[[locID.ind]]) # problem
              }
              check.ind <- random.points.to.eval %in% row.names(glmssn$sampinfo$X)
              random.points.to.eval <- random.points.to.eval[check.ind]
              d.list.to.eval[[j]] <- c(random.points.to.eval, extra.arguments$fixed)
              d.list[[j]] <- c(random.points, fixed.points)
            } else {
              d.list.to.eval[[j]] <- d.list[[j]] <- c(random.points, extra.arguments$fixed)
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
      
      # Print progress to user
      
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
  
  # Update SSN object
  
  final.points <<- unique(final.points)
  # suppressWarnings(subsetSSN(ssn = ssn, filename = new.ssn.path, pid %in% final.points))
  # preds <- NULL
  # if(length(ssn@predpoints@SSNPoints) > 0){
  #   preds <- "preds"
  # }
  # ssn.new <- importSSN(new.ssn.path, preds)
  # createDistMat(ssn.new, preds, TRUE, TRUE)
  # rm(final.points, pos = 1) # For some reason, the subsetSSN doesn't work unless final.points is in the global scope. Removing it here.
  # 
  # return updated ssn
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
