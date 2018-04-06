#' Get row or element numbers for matches between two vectors
#' 
#'@description
#'
#'Get the positions in a vector corresponding to a match between vectors.
#'
#'@usage 
#'
#'\code{matchIndices} is a function that returns row or element numbers where a set of values from one object occur in another. 
#'\code{matchIndices(match.from, match.to)}. 
#'
#'Note, this function is only intended to be used internally by \code{generateSites}.
#'
#'@param match.from a numeric or character vector 
#'@param match.to a numeric or character vector
#'@return A numeric vector containing the row indices of every match, in the order that they appear in match.from
#' 
#'@export
matchIndices <- function(match.from, match.to){
  n1 <- length(match.from)
  n2 <- length(match.to)
  if(n1 != n2){
    stop("Vectors must be of the same length.")
  }
  match.from <- as.character(match.from)
  match.to <- as.character(match.to)
  match.made <- rep(NA, n1)
  for(i in 1:n1){
    match.now <- match.from[i]
    match.made[i] <- which(match.to == match.now)
  }
  if(sum(is.na(match.made)) > 0){
    stop("Matching could not be complete.")
  }
  return(match.made)
}
