#' A utility function for parameter estimation with known covariance parameters (sequential)
#' 
#' @description
#'   
#' The function \code{sequentialDOptimality} is a utility function that can be used with \code{\link{optimiseSSNDesign}}. It is a utility function that minimises the determinant of the variance-covariance matrix of the fixed effects for models that are fit over a set of design points. 
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param glmssn An model object of class glmssn.
#' @param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#' @param prior.parameters A list of random functions that are parameterised in terms of n.draws.
#' @param n.draws A numeric scalar for the number of Monte Carlo draws to use when approximating the utility. 
#' @param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{optimiseSSNDesign}}.
#' @return A numeric scalar.
#' 
#' @details
#'
#' Note that this function operates differently to \code{DOptimality}. The functions \code{DOptimality} and \code{EDOptimality} assume there are no sites which have already been incorporated into a design. They compute the variance-covariance matrix on the fixed effects (Som et al., 2014) for each set of simulated or esimated covariance parameters, respectively. In this sequential form, the observed variance-covariance matrix is extracted for the sites which are fixed in the design. Then the sites which are not fixed are used to estimate the expected variance-covariance matrix. The observed and expected matrices are added, before being reduced to a scalar value to serve as the utility.
#' 
#' @export
sequentialDOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Extract out the variance-covariance matrix for the fixed effects from the glmssn object
  old.covbi <- glmssn$estimates$covbi
  
  # Retrieve relevant design matirx
  ind <- row.names(extra.arguments$obs.X) %in% design.points
  X <- extra.arguments$obs.X[ind,]
  Xt <- t(X)
  
  # Retrieve coordinates of sampling points
  ind <- row.names(extra.arguments$obs.C) %in% design.points
  cds <- extra.arguments$obs.C[ind,]
  
  # Cut down distance matrices, etc.
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  n.zero <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
  ## Get other model parameters
  td <- glmssn$args$useTailDownWeight
  cm <- glmssn$args$CorModels
  un <- glmssn$args$use.nugget
  ua <- glmssn$args$use.anisotropy
  re <- glmssn$sampInfo$REs
  
  ## Simulate covparms from priors
  cvp.cols <- length(glmssn$estimates$theta)
  cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  for(i in 1:cvp.cols){
    cvp[, i] <- prior.parameters[[i]](n.draws)
  }
  
  ## Perform MC simulations
  D <- vector("numeric", n.draws)
  
  for(i in 1:n.draws){
    
    theta.i <- cvp[i, ]
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
  #  covb <- solve(covbi)
    D[i] <- log(det(covbi + old.covbi))
    
  }
  
  return(mean(D, na.rm = TRUE))
  
}
