#' Search all possible designs to find one which maximises the utility
#' 
#'@description
#'
#'This function takes all possible combinations of n design points and evaluates a utility function for each combination.
#' 
#'@usage 
#'
#'\code{searchAllDesigns(ssn, glmssn, afv.column, n.points, utility.function, prior.parameters, n.draws = 500, extra.arguments = NULL)}
#' 
#'@param ssn an object of class SpatialStreamNetwork
#'@param glmssn an object of class glmssn
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
#'This is very computationally expensive. Only use this function for small examples or to benchmark the performance of the greedy exchange algorithm.
#'
#'@export
searchAllDesigns <- function(ssn, n.points, model, utility.function, prior.parameters, n.draws = 500, extra.arguments = NULL){
  
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
  if(!is.list(extra.arguments)){
    extra.arguments <- list()
  }
  
  # Extract distance and connectivity matrices for each network
  
  # Extract important matrices for each network
  
  dist.junc.obs <- getStreamDistMatInOrder(ssn)
  if(length(ssn@predpoints@SSNPoints) > 0){
    dist.junc.prd <- getStreamDistMatInOrder(ssn, "preds")
    dist.junc.pxo <- getStreamDistMatInOrder.predsxobs(ssn)
    network.each.pred <- ssn@predpoints@SSNPoints[[1]]@point.data$netID
  }
  
  # Loop for each network
  
  w.in.network.maxes <- vector("list", n)
  for(net in 1:n){
    
    # Find network specific data
    
    n.final.points.this.network <- n.points[net]
    # dist <- dc.matrices$Distance[[net]]
    # conn <- dc.matrices$Connectivity[[net]]
    ind <- ssn@obspoints@SSNPoints[[1]]@point.data$netID == net
    points.this.network <- ssn@obspoints@SSNPoints[[1]]@point.data$pid[ind]
    designs.this.network <- t(combn(points.this.network, n.final.points.this.network))
    
    # Get network's distance matrices
    
    dist.junc.obs.net <- dist.junc.obs[[net]]
    
    if(length(ssn@predpoints@SSNPoints) > 0){
      dist.junc.prd.net <- dist.junc.prd[[net]]
      dist.junc.pxo.a.net <- dist.junc.pxo[[ 2 * net - 1 ]]
      dist.junc.pxo.b.net <- dist.junc.pxo[[ 2 * net ]]
    }
    
    extra.arguments$Matrices.Obs <- getImportantMatrices.obs(
      dist.junc.obs.net, 
      afv = ssn@obspoints@SSNPoints[[1]]@point.data[ind,afv.column]
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
        afv.obs = ssn@obspoints@SSNPoints[[1]]@point.data[ind,afv.column],
        afv.prd = ssn@predpoints@SSNPoints[[1]]@point.data[indp,afv.column]
      )
    }
    
    # Evaluate the utility for all combinations
    
    design.data <- ssn@obspoints@SSNPoints[[1]]@point.data
    n.iter <- nrow(designs.this.network)
    U <- vector("numeric", n.iter)
    for(i in 1:n.iter){
      
      # Extract out design and data associated with design
      design.i <- designs.this.network[i, ]
      #design.data.i <- design.data[design.data$locID %in% design.i, ]
      # ind <- row.names(dist) %in% design.i
      #dist.i <- dist[ind, ind]
      #conn.i <- conn[ind, ind]
        
        U[i] <- utility.function(
          ssn,
          glmssn,
          design.i,
          prior.parameters, 
          n.draws,
          extra.arguments
        )
      
    }
    
    # Find and record design that maximises utility
    
    indm <- U == max(U)
    w.in.network.maxes[[net]] <- designs.this.network[indm, ]
    
  }
  
  ## create new SSN containing only the selected points
  final.points <- unlist(w.in.network.maxes)
  
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
