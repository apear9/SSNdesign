identifyObsxPredsABMatrices <- function(
  dist.matrices
){
  
  # Get names
  
  dist.names <- names(dist.matrices)
  
  # Find names corresponding to left-and-right obsxpreds matrices
  
  ind1 <- grep(".a", dist.names)
  ind2 <- grep(".b", dist.names)
  inds <- union(ind1, ind2)
  
  # Spit out indices
  
  return(inds)
  
}