#' A function to crush down an edgelist in to a vector of paired vertices
#' 
#' @param edgelist An n x 2 matrix, as produced by \code{get.edgelist()}
#' @return A vector. 
collapseAndReorderEdgelist <- function(edgelist){
  squashed_edgelist <- matrix(edgelist, ncol = 1)
  nrows <- nrow(edgelist)
  ordering <- c(2 * 1:nrows - 1, 2 * 1:nrows)
  squashed_edgelist <- cbind(squashed_edgelist, ordering)
  squashed_edgelist <- squashed_edgelist[order(squashed_edgelist[,"ordering"]), -2]
  return(squashed_edgelist)
}