#' A function to solve adaptive / sequential design problems on SpatialStreamNetworks
#'@description 
#'
#'Finds an optimal design consisting of a specified number of points on a SpatialStreamNetwork, given a set of sites that must remain fixed. 
#' 
#'@usage 
#' 
#'\code{doAdaptiveDesign(ssn, glmssn, fixed.points, afv.column, n.points, utility.function, prior.parameters, n.draws = 500, extra.arguments = NULL)}
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
#'@param extra.arguments a list of any extra parameters which can be used to control the behaviour of this function or the utility function
#'@return An object of class SpatialStreamNetwork. The SSNPoints for the obspoints slot will be updated to reflect the selected design. 
#'
#'@details 
#'
#'This function implements the greedy exchange algorithm per Falk et al. (2014) to quickly find an optimal design from a set of sites. The algorithm does not always discover the absolute best design; however, even when this occurs, the algorithm will quickly discover a reasonably effective design. 
#'  
#'@export
doAdaptiveDesign <- function(
  ssn,
  new.ssn.path,
  glmssn,
  fixed.points,
  afv.column, 
  n.points, 
  utility.function, 
  prior.parameters, 
  n.draws = 500, 
  extra.arguments = NULL){
  
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
  ### PARAMETER TO CONTROL WHETHER DESIGNS ARE FOUND BY NETWORK OR ACROSS NETWORKS
  do.separately <- TRUE
  if(length(n.points) == 1){
    do.separately <- FALSE
  }
  
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
  
  glmssn$sampinfo$X <- model.matrix(glmssn$args$formula, data = ssn@obspoints@SSNPoints[[1]]@point.data)
  
  # Extract important matrices for each network
  
  dist.junc.obs <- getStreamDistMatInOrder(ssn)
  if(length(ssn@predpoints@SSNPoints) > 0){
    dist.junc.prd <- getStreamDistMatInOrder(ssn, "preds")
    dist.junc.pxo <- getStreamDistMatInOrder.predsxobs(ssn)
  }
  
  # Initialise for loop
  
  all.points <- ssn@obspoints@SSNPoints[[1]]@point.data$pid
  
  #### SPLIT HERE TO EITHER RUN EACH NETWORK SEPARATELY OR RUN ALL SITES TOGETHER
  
  if(do.separately){
    
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
      
      # vector and list for holding designs in each of K iterations
      
      U.all <- c()
      design.all <- list()
      for(k in 1:K){
        
        # Select random subsample of points
        random.sample <- sample(1:n.points.this.network, n.final.this.network, FALSE)
        random.points <- points.this.network[random.sample]
        random.points <- as.character(random.points)
        random.and.fixed.points <- c(random.points, as.character(fixed.points.this.network))
        
        # Evaluate utility once to initialise
        U <- utility.function(
          ssn,
          glmssn,
          random.and.fixed.points,
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
            
            for(j in 1:n.j){
              
              random.points[i] <- remaining.values[j]
              random.and.fixed.points <- c(random.points, as.character(fixed.points.this.network))
              designs <- append(designs, list(random.and.fixed.points))
              U_ <- utility.function(
                ssn,
                glmssn,
                random.and.fixed.points,
                prior.parameters, 
                n.draws,
                extra.arguments
              )
              
              U.ij <- append(U.ij, U_)
              
            }
            
          }
          
          design.previous <- design.now
          if(max(U.ij) >= U){
            U <- max(U.ij)
            ind <- U.ij == max(U.ij)
            ind <- (1:length(U.ij))[ind]
            #U.ij <- U.ij[-ind]
            design.now <- designs[[ind]]
            random.points <- design.now
            
          }
          
          cond <- !all(design.now %in% design.previous)
          #cnt = cnt + 1
        }
        # print(cnt)
        U.all <- append(U.all, U)
        design.all <- append(design.all, list(design.now))
        
      }
      
      ind <- (1:length(U.all))[U.all == max(U.all)][1]
      final.points[[net]] <- design.all[[ind]]
    }
    
    ## create new SSN containing only the selected points
    
    final.points <- unlist(final.points)
    
  } else {
    
    extra.arguments$fixed <- fixed.points
    
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
      
      # Select random subsample of points
      random.sample <- sample(1:length(all.points), n.points, FALSE)
      random.points <- all.points[random.sample]
      random.and.fixed.points <- c(random.points, fixed.points)
      random.and.fixed.points <- as.character(random.and.fixed.points)
      
      # Evaluate utility once to initialise
      
      U <- utility.function(
        ssn,
        glmssn,
        random.and.fixed.points,
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
          
          for(j in 1:n.j){
            
            random.points[i] <- remaining.values[j]
            random.and.fixed.points <- c(random.points, as.character(fixed.points))
            designs <- append(designs, list(random.and.fixed.points))
            U_ <- utility.function(
              ssn,
              glmssn,
              random.and.fixed.points,
              prior.parameters, 
              n.draws,
              extra.arguments
            )
            
            U.ij <- append(U.ij, U_)
            
          }
          
        }
        
        design.previous <- design.now
        if(max(U.ij) >= U){
          U <- max(U.ij)
          ind <- U.ij == max(U.ij)
          ind <- (1:length(U.ij))[ind]
          design.now <- designs[[ind]]
          random.points <- design.now
          
        }
        
        cond <- !all(design.now %in% design.previous)
        
      }
      
      U.all <- append(U.all, U)
      design.all <- append(design.all, list(design.now))
      
    }
    
    ind.max <- U.all == max(U.all)
    final.points <- unlist(design.all[ind.max])
    
  }
  
  # # update point coords
  # ind.point.coords <- attributes(ssn@obspoints@SSNPoints[[1]]@point.coords)$dimnames[[1]] %in% final.points
  # ssn@obspoints@SSNPoints[[1]]@point.coords <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind.point.coords, ]
  # 
  # # update network point coords
  # ind.network.point.coords <- row.names(ssn@obspoints@SSNPoints[[1]]@network.point.coords) %in% final.points
  # ssn@obspoints@SSNPoints[[1]]@network.point.coords <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.network.point.coords, ]
  # 
  # # update point data
  # ind.point.data <- row.names(ssn@obspoints@SSNPoints[[1]]@point.data) %in% final.points
  # ssn@obspoints@SSNPoints[[1]]@point.data <- ssn@obspoints@SSNPoints[[1]]@point.data[ind.point.data, ]
  # 
  # # update bbox
  # new.bbox <- sp::bbox(ssn@obspoints@SSNPoints[[1]]@point.coords)
  # ssn@obspoints@SSNPoints[[1]]@points.bbox <- new.bbox
  final.points <<- final.points
  suppressWarnings(subsetSSN(ssn = ssn, filename = new.ssn.path, pid %in% final.points))
  preds <- NULL
  if(length(ssn@predpoints@SSNPoints) > 0){
    preds <- "preds"
  }
  ssn.new <- importSSN(new.ssn.path, preds)
  createDistMat(ssn.new, preds, TRUE, TRUE)
  rm(final.points, pos = 1) # For some reason, the subsetSSN doesn't work unless final.points is in the global scope. Removing it here.
  
  # return updated ssn
  return(ssn.new)
  
}
