#' Get the matrices required to calculate the covariance matrix on the data
#' 
#'@description
#'
#'Calculates the distance and related matrices between observed or prediction sites on a SpatialStreamNetwork.
#'
#'@usage
#'
#'\code{getImportantMatrices(d.junc, afv = NULL)}
#'    
#'@param d.junc An asymmetric, square distance matrix on the design points. 
#'@param afv The additive function values for each of the design points, in order of their pid.
#'@return A list of 5 matrices in the order dist.hydro, a.mat, b.mat, conn.mat, and w.mat.
#'
#'@details
#'
#'This function is called internally by \code{\link{findOptimalDesign}} and \code{\link{doAdaptiveDesign}}. This is done to create the distance and related matrices required to calculate the covariances between observed sites, which is needed for all the utility functions. Note that, although the name of the function suggests this should only work for observed sites, it also works for prediction sites. Therefore, it is also called by \code{\link{findOptimalDesign}} and \code{\link{doAdaptiveDesign}} to create the distance and related matrices among the prediction sites.
#'
#'@export 
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
