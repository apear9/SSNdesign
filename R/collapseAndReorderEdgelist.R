#' A function to crush down an edgelist in to a vector of paired vertices
#'
#'@description
#'
#'An edgelist is an n x 2 matrix which details connections between edges in a network; from the first column to the second column. This function reduces the matrix to a vector which is more convenient to deal with inside the function \code{generateSites}.
#'
#'@usage
#'
#'\code{collapseAndReorderEdgelist(edgelist)}
#'         
#'@param edgelist An n x 2 matrix, as produced by the function \code{get.edgelist()}
#'@return A numeric vector.
#'
#'@details
#'
#'This is a function called inside \code{generateSites} to correctly place simulated sites on a SpatialStreamNetwork. It is not intended for use outside of this context, though it may prove useful in other applications.
#'   
#'@export 
collapseAndReorderEdgelist <- function(edgelist){
  squashed_edgelist <- matrix(edgelist, ncol = 1)
  nrows <- nrow(edgelist)
  ordering <- c(2 * 1:nrows - 1, 2 * 1:nrows)
  squashed_edgelist <- cbind(squashed_edgelist, ordering)
  squashed_edgelist <- squashed_edgelist[order(squashed_edgelist[,"ordering"]), -2]
  return(squashed_edgelist)
}