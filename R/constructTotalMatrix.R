constructTotalMatrix <- function(list.matrices, pxo = FALSE){
  
  # Input checking
  if(!is.list(list.matrices)){
    stop("The argument list.matrices must be a list of matrices.")
  }
  if(pxo){
    keep.ind <- grep("a", names(list.matrices))
    list.matrices <- list.matrices[keep.ind]
  }
  
  # Find size of list as number of networks
  n <- length(list.matrices) 
  
  # Find the dimensions of the output matrix
  rows <- unlist(lapply(list.matrices, nrow))
  cols <- unlist(lapply(list.matrices, ncol))
  nrows <- sum(rows)
  ncols <- sum(cols)
  
  # Pre-generate the empty matrix
  total.matrix <- net.zero <- matrix(0, nrows, ncols)
  
  # Find the cumulative sum over the rows and columns as starting indices
  crows <- cumsum(rows)
  ccols <- cumsum(cols)
  
  # Define last.row and last.col as anchors for the filling in
  
  last.row <- 1
  last.col <- 1
  
  # Loop through the matrix for each network and fill in the matrix as required.
  
  for(i in 1:n){
    
    total.matrix[last.row:crows[i], last.col:ccols[i]] <- list.matrices[[i]]
    net.zero[last.row:crows[i], last.col:ccols[i]] <- 1
    
    last.row <- 1 + crows[i]
    last.col <- 1 + ccols[i]
    
  }
  
  # Attach names to matrix rows and columns
  rnames <- unlist(lapply(list.matrices, rownames))
  cnames <- unlist(lapply(list.matrices, colnames))
  rownames(total.matrix) <- rnames
  colnames(total.matrix) <- cnames
  rownames(net.zero) <- rnames
  colnames(net.zero) <- cnames
  
  # Return the outputs
  return(list(d.junc = total.matrix, net.zero = net.zero))
  
}
