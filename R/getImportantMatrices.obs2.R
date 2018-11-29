getImportantMatrices.obs2 <- function(d.junc, afv = NULL){
  
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
