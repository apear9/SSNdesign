getStreamDistMatInOrder.predsxobs <- function(x){
  
  # Get stream distance matrices
  
  unordered.distmats <- getStreamDistMat(x, "preds")
  
  # If for preds, identify A, B matrices and remove
  
  ind <- identifyObsxPredsABMatrices(unordered.distmats)
    
  unordered.distmats <- unordered.distmats[ind]
    
  # Names are in alphabetical order
  
  unordered.names <- names(unordered.distmats)
  
  # Strip away the parts of the names that include dist.net
  
  unordered.names <- gsub("dist.net", "", unordered.names)
  
  # Strip away the parts that include .a and .b
  
  unordered.names <- gsub(".[ab]{1}", "", unordered.names)
  
  # Convert to numbers
  
  unordered.names <- as.numeric(unordered.names)
  
  # Find correct sorting order
  
  sort.order <- order(unordered.names)
  
  # Sort list
  
  ordered.distmats <- unordered.distmats[sort.order]
  
  # Return result
  
  return(ordered.distmats)
  
}
