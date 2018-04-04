#' Implement optimal designs on a SpatialStreamNetwork
#'
#' \code{findOptimalDesign()} 
#' 
#'  @param ssn an object of class SpatialStreamNetwork
#'  @param glmssn an object of class glmssn
#'  @param afv.column the name of the column in the SpatialStreamNetwork object that contains the additive function values
#'  @param n.points the number of points to be included in the final design. This can be a single number, in which case all networks will have the same number of points. This can also be a vector with the same length as the number of networks in ssn. 
#'  @param utility.function a function with the signature Utility Function. Built-in functions are Doptimality, FisherInformationMatrix, ... 
#'  @param prior.parameters a function to act as a prior for covariance parameter values
#'  @param n.draws a numeric value, being the number of Monte Carlo draws to take when evaluating potential designs
#'  @param extra.arguments a list of any extra parameters (see below) which can be used to control the behaviour of this function or the utility function
#'  @return An object of class SpatialStreamNetwork with its observed SSNPoints slot updated to reflect the optimal design selected under the given utility function.
#'  
#'  @export
findOptimalDesign <- function(
  ssn, 
  glmssn, 
  afv.column, 
  n.points, 
  utility.function, 
  prior.parameters, 
  n.draws = 500, 
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
  if(!is.numeric(n.points) | n.points < 1){
    stop("The argument n.points must have a positive integer value.")
  }
  if(length(n.points) == 1){
    n.points <- rep(n.points, n)
  }
  if(!is.numeric(n.draws) | n.draws < 1){
    stop("The argument n.draws must have a positive integer value.")
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
  
  # Extract important matrices for each network
  
  dist.junc.obs <- getStreamDistMatInOrder(ssn)
  if(length(ssn@predpoints@SSNPoints) > 0){
    dist.junc.prd <- getStreamDistMatInOrder(ssn, "preds")
    dist.junc.pxo <- getStreamDistMatInOrder.predsxobs(ssn)
  }
  
  # Initialise for loop
  
  all.points <- ssn@obspoints@SSNPoints[[1]]@point.data$pid
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
    n.points.this.network <- length(points.this.network)
    
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
      extra.arguments$Matrices.pxo <- getImportantMatrices.pxo(
        dist.junc.pxo.a.net,
        dist.junc.pxo.b.net,
        afv.obs = ssn@obspoints@SSNPoints[[1]]@point.data[ind1,afv.column],
        afv.prd = ssn@predpoints@SSNPoints[[1]]@point.data[indp,afv.column]
      )
    }
    
    # vector and list for holding designs in each of K iterations
    
    U.all <- c()
    design.all <- list()
    for(k in 1:K){
      
      # Select random subsample of points
      
      ## BREAK POINT -- does not work if pred sites present
      
      # ind1 <- ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID == net
      # points.this.network <- ssn@obspoints@SSNPoints[[1]]@point.data$pid[ind1]
      # points.this.network <- as.numeric(points.this.network)
      # n.points.this.network <- length(points.this.network)
      random.sample <- sample(1:n.points.this.network, n.final.this.network, FALSE)
      random.points <- points.this.network[random.sample]
      random.points <- as.character(random.points)
      # print(random.points)
      
      # Evaluate utility once to initialise
      # OBSOLETE CODE -- need to define distance matrices differently
      # dist.tot <- dc.matrices$Distance[[net]]
      # conn.tot <- dc.matrices$Connectivity[[net]]
      # ind <- row.names(dist.tot) %in% random.points
      # dist.ij <- dist.tot[ind, ind]
      # conn.ij <- conn.tot[ind, ind]
      # design.data <- ssn@obspoints@SSNPoints[[1]]@point.data
      # design.data <- design.data[order(design.data$locID), ]
      # design.data.ij <- design.data[design.data$locID %in% random.points, ]
      # 
      U <- utility.function(
        ssn,
        glmssn,
        random.points,
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
      
      #cnt = 0
      
      while(cond){
        
        for(i in 1:n.final.this.network){
          
          last.value <- random.points
          ind.i <- !(points.this.network %in% last.value)
          remaining.values <- points.this.network[ind.i]
          remaining.values <- as.character(remaining.values)
          n.j <- length(remaining.values)
          
          for(j in 1:n.j){
            
            random.points[i] <- remaining.values[j]
            designs <- append(designs, list(random.points))
            # ind <- row.names(dist.tot) %in% random.points
            # dist.ij <- dist.tot[ind, ind]
            # conn.ij <- conn.tot[ind, ind]
            # design.data.ij <- design.data[design.data$locID %in% random.points, ]
            # 
            U_ <- utility.function(
              ssn,
              glmssn,
              random.points,
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
  
  # update point coords
  ind.point.coords <- attributes(ssn@obspoints@SSNPoints[[1]]@point.coords)$dimnames[[1]] %in% final.points
  ssn@obspoints@SSNPoints[[1]]@point.coords <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind.point.coords, ]
  
  # update network point coords
  ind.network.point.coords <- row.names(ssn@obspoints@SSNPoints[[1]]@network.point.coords) %in% final.points
  ssn@obspoints@SSNPoints[[1]]@network.point.coords <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.network.point.coords, ]
  
  # update point data
  ind.point.data <- row.names(ssn@obspoints@SSNPoints[[1]]@point.data) %in% final.points
  ssn@obspoints@SSNPoints[[1]]@point.data <- ssn@obspoints@SSNPoints[[1]]@point.data[ind.point.data, ]
  
  # update bbox
  new.bbox <- sp::bbox(ssn@obspoints@SSNPoints[[1]]@point.coords)
  ssn@obspoints@SSNPoints[[1]]@points.bbox <- new.bbox
  
  # return updated ssn
  return(ssn)
  
}
