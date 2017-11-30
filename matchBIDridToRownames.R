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
