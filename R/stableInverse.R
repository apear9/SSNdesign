#' A function to compute inverses by eigenvalue decomposition, avoiding computational singularities
#' 
#' \code{stableInverse()}
#' 
#' @param mat A square matrix object whose inverse is desired.
#' @param threshold Eigenvalues below this threshold will be dropped and an approximate inverse will be returned. 
stableInverse <- function(mat, threshold = 1e-10){
  
  eig <- eigen(mat)
  evl <- eig$values
  evc <- eig$vectors
  evct <- t(evc)
  ind <- evl > threshold
  evl <- diag(1/evl[ind])
  evc <- evc[,ind]
  evct <- evct[ind,]
  return(evc %*% evl %*% evct)
  
}
