#' A utility function for fixed effects parameter estimation with known covariance parameters
#' 
#'@description
#'
#'\code{DOoptimality} is a utility function that can be used with \code{\link{findOptimalDesign}}. It is a utility function that minimises the determinant of the variance-covariance matrix of the fixed effects over a set of designs. 
#' 
#'@usage
#'
#'\code{DOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
#' 
#'@param ssn An object of class SpatialStreamNetwork
#'@param glmssn An model object of class glmssn.
#'@param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#'@param prior.parameters A list of random functions that are parameterised in terms of n.draws.
#'@param n.draws A numeric scalar for the number of Monte Carlo draws to use when approximating the utility. 
#'@param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{findOptimalDesign}} and \code{\link{doAdaptiveDesign}}.
#'@return A numeric scalar.
#' 
#'@details
#'
#'This utility function is similar to \code{\link{EDOptimality}}, except that it assumes the covariance parameters are known. The definition of this utility function follows that in Som et al. (2014), with a minor modification to allow the incorporation of prior information in the form of Monte Carlo simulations.
#'
#' @export
DOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){

  ## Get data for design points
  # Rely on ordering of output of glmssn to retrieve correct rows in the design matrix
  ind <- row.names(glmssn$sampinfo$X) %in% design.points
  X <- glmssn$sampinfo$X[ind, ]
  Xt <- t(X)
  
  ind.cds <- row.names(ssn@obspoints@SSNPoints[[1]]@point.coords) %in% design.points
  cds <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind.cds, ]
  colnames(cds) <- c("x", "y")
  
  ## Get distance matrices, etc.
  
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  n.zero <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
  
  ## Get simulated covariance parameters
  
  cvp.cols <- length(glmssn$estimates$theta)
  cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  
  ## Get other model parameters
  
  td <- glmssn$args$useTailDownWeight
  cm <- glmssn$args$CorModels
  un <- glmssn$args$use.nugget
  ua <- glmssn$args$use.anisotropy
  re <- glmssn$sampInfo$REs
  
  ## Simulate covparms from priors
  
  for(i in 1:cvp.cols){
    
    cvp[, i] <- prior.parameters[[i]](n.draws)
    
  }
  
  ## Perform MC simulations
  
  D <- vector("numeric", n.draws)
  
  for(i in 1:n.draws){
    
    theta.i <- cvp[i, ]
    #print(n.zero)
    V <- SSN:::makeCovMat(
      theta.i, 
      mat$d, 
      mat$a, 
      mat$b, 
      mat$w, 
      n.zero, 
      cds[, "x"], 
      cds[, "y"], 
      cds[, "x"], 
      cds[, "y"], 
      td,
      cm, 
      un,
      ua,
      re
    )
    
    covbi <- t(X) %*% solve(V) %*% X
    covb <- solve(covbi)
    D[i] <- -log(det(covb))
    
  }
  
  return(mean(D, na.rm = TRUE))
  
}
