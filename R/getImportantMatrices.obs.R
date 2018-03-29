#' Get the matrices required to calculate the covariance matrix on the data
#' 
#' \code{getImportantMatrices()}
#' 
#' @param d.junc A non-symmetric distance matrix on the design points. See... for more details. 
#' @param afv The additive function values for each of the design points, in order of their pid.
#' @return A list of 5 matrices in the order dist.hydro, a.mat, b.mat, conn.mat, and w.mat.
#' 
getImportantMatrices.obs <- function(d.junc, afv = NULL){
  
  # Input checking 
  if(!is.matrix(d.junc)){
    stop("The input d.junc must a matrix.")
  }
  
  # Extract hydrologic, a, and b-distance matrices
  
  d.hydro <- d.junc + t(d.junc)
  a.mat <- pmax(d.junc, t(d.junc))
  b.mat <- pmin(d.junc, t(d.junc))
  
  # From b-distance matrix, derive connectivity matrix
  
  c.mat <- 1 - 1 * (b.mat > 0)


  # Create weights matrix
  n.mat <- nrow(d.hydro)
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
  w.mat <- c.mat * sqrt(
    pmin(
      outer(afv,rep(1, times = n.mat)), 
      t(outer(afv,rep(1, times = n.mat)))
    ) /
      pmax(
        outer(afv,rep(1, times = n.mat)),
        t(outer(afv,rep(1, times = n.mat)))
      )
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
