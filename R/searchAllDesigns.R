#' Search all possible designs to find one which maximises the utility
#' 
#' \code{searchAllDesigns()}
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param n.points n.points the number of points to be included in the final design. This can be a single number, in which case all networks will have the same number of points. This can also be a vector with the same length as the number of networks in ssn. 
#' @param model either a formula object or a linear model object (of classes lm, glm, etc.) from which a model formula can be extracted.
#' @param utility.function a function with the signature Utility Function. Built-in functions are Doptimality, FisherInformationMatrix, ... 
#' @param prior.parameters a function to act as a prior for covariance parameter values
#' @param n.draws a numeric value, being the number of Monte Carlo draws to take when evaluating potential designs
#' @param extra.arguments a list containing any extra arguments that must be passed to the utility function
#' @return An object of class SpatialStreamNetwork with its SSNPoints slot updated to reflect the optimal design selected under the given utility function.
#'
#'
searchAllDesigns <- function(ssn, n.points, model, utility.function, prior.parameters, n.draws = 500, extra.arguments){
  
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
  
  
  # Extract distance and connectivity matrices for each network
  
  dc.matrices <- retrieveDistAndConnMat(ssn)
  
  # Loop for each network
  
  w.in.network.maxes <- vector("list", n)
  for(net in 1:n){
    
    # Find network specific data
    
    n.final.points.this.network <- n.points[net]
    dist <- dc.matrices$Distance[[net]]
    conn <- dc.matrices$Connectivity[[net]]
    ind <- ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID == net
    points.this.network <- row.names(ssn@obspoints@SSNPoints[[1]]@network.point.coords)[ind]
    designs.this.network <- t(combn(points.this.network, n.final.points.this.network))
    
    # Evaluate the utility for all combinations
    
    design.data <- ssn@obspoints@SSNPoints[[1]]@point.data
    n.iter <- nrow(designs.this.network)
    U <- vector("numeric", n.iter)
    for(i in 1:n.iter){
      
      # Extract out design and data associated with design
      design.i <- designs.this.network[i, ]
      design.data.i <- design.data[design.data$locID %in% design.i, ]
      ind <- row.names(dist) %in% design.i
      dist.i <- dist[ind, ind]
      conn.i <- conn[ind, ind]
        
        U[i] <- utility.function(
          formula = model, 
          design = design.data.i, 
          important.matrices = list(Distance = dist.i, Connectivity = conn.i), 
          prior.parameters = prior.parameters, 
          n.draws = n.draws,
          extra.arguments = extra.arguments
        )
      
    }
    
    # Find and record design that maximises utility
    
    ind <- U == max(U)
    w.in.network.maxes[[net]] <- designs.this.network[ind, ]
    
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
