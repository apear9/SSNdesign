stableInverse <- function(mat, threshold = 1e-10){
  
  # Compute eigenvalue decomposition
  eig <- eigen(mat)
  evl <- eig$values # Eigenvalues
  evc <- eig$vectors # Eigenvectors in the matrix V (column space)
  evct <- t(evc) # Eigenvectors in the matrix V^T (row space)
  # Remove eigenvalues below the threshold
  ind <- evl > threshold 
  # Invert the diagonal matrix of the remaining eigenvalues
  evl <- diag(1/evl[ind])
  # Remove the corresponding eigenvectors
  evc <- evc[,ind]
  evct <- evct[ind,]
  # The result will be an approximation to the inverse (non-unique) with the same dimensions as the original matrix
  return(evc %*% evl %*% evct)
  
}
