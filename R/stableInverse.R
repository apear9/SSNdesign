#' A function to compute inverses by eigenvalue decomposition, avoiding computational singularities
#' 
#'@description
#'
#'\code{stableInverse} is a function designed to approximate the inverse of a square matrix which is computationally singular. It does this by using the eigenvalue decomposition of the matrix.
#' 
#'@param mat A square matrix object whose inverse is desired.
#'@param threshold Eigenvalues below this threshold will be dropped and an approximate inverse will be returned. 
#'@return The inverse of the input matrix.
#'
#'@details
#'
#'In inverting the Fisher information matrix, e.g. as required by the \code{\link{CPOptimality}} utility function, it often happens that one or more of the eigenvalues are very close to zero. This results in a computational singularity; i.e. the matrix cannot be inverted. This may occur in spite of the Fisher information matrix being analytically symmetric positive definite. In this case, the eigenvalue decomposition of the matrix is computed. Any eigenvalues below a threshold (by default, 1e-10) will be omitted and so will the corresponding eigenvectors in the diagonalised matrix equation VUV' where V is the matrix containing the eigenvectors as columns and U is the diagonal matrix of the eigenvalues. This matrix equation will still return a matrix of the same dimensions as the original matrix. This resulting matrix will be approximately equal to the inverse of the original matrix.
#'   
#'@export
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
