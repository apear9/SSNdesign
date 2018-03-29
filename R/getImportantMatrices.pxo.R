#' Get the matrices required to calculate the covariance matrix on the data
#' 
#' \code{getImportantMatrices()}
#' 
#' @param d.junc A non-symmetric distance matrix on the design points. See... for more details. 
#' @param afv The additive function values for each of the design points, in order of their pid.
#' @return A list of 5 matrices in the order dist.hydro, a.mat, b.mat, conn.mat, and w.mat.
#' 
getImportantMatrices.pxo <- function(d.junc.a, d.junc.b, afv.obs = NULL, afv.prd = NULL){
  
  # Input checking 
  if(!is.matrix(d.junc.a) | !is.matrix(d.junc.b)){
    stop("The inputs d.junc.a and d.junc.b must be a matrix.")
  }
  
  # Extract hydrologic, a, and b-distance matrices
  
  d.hydro <- as.matrix(d.junc.a + t(d.junc.b))
  a.mat <- pmax(d.junc.a, t(d.junc.b))
  b.mat <- pmin(d.junc.a, t(d.junc.b))
  
  # From b-distance matrix, derive connectivity matrix
  
  c.mat <- 1 - 1 * (b.mat > 0)
  
  # Create weights matrix
  n.prd <- ncol(d.hydro)
  n.obs <- nrow(d.hydro)
  # print(c.mat)
  # print(n.mat == ncol(d.hydro))
  # print(sqrt(
  #   pmin(
  #     outer(afv,rep(1, times = n.mat)), 
  #     t(outer(afv,rep(1, times = n.mat)))
  #   ) /
  #     pmax(
  #       outer(afv,rep(1, times = n.mat)),
  #       t(outer(afv,rep(1, times = n.mat)))
  #     )
  # ))
  #print(c.mat)
  w.mat <- c.mat * sqrt(
    pmin( outer( afv.obs, rep(1, times = n.prd) ), t( outer( afv.prd, rep(1, times = n.obs) ) ) ) 
    /
    pmax( outer( afv.obs, rep(1, times = n.prd) ), t( outer( afv.prd, rep(1, times = n.obs) ) ) )
  )
  
  # Return all matrices
  
  return(
    list(
      d = d.hydro,
      a = a.mat,
      b = b.mat,
      c = c.mat,
      w = w.mat
    )
  )
  
}
