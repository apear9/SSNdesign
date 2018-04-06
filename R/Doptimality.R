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
  
  # OBSOLETE CODE
  # ## Get distance matrices, etc.
  # 
  # distance.matrix <- important.matrices$Distance
  # connectivity.matrix <- important.matrices$Connectivity
  # weights.vector <- rep(1, nrow(connectivity.matrix))
  # conn.weight.matrix <- connectivity.matrix 
  # 
  # ## Get other arguments
  # 
  # cvp <- vector("list", length(prior.parameters)) 
  # n <- nrow(distance.matrix)
  # X <- model.matrix(object = formula, data = design)
  # Xt <- t(X)
  
  ## Get data for design points
  
  # print(design.points)
  
  # Rely on ordering of output of glmssn to retrieve correct rows in the design matrix
  ind <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  
  # print(ind)
  
  X <- glmssn$sampinfo$X[ind, ]
  
  # print(X)
  
  Xt <- t(X)
  cds <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind, ]
  # print(class(cds))
  # 
  colnames(cds) <- c("x", "y")
  
  ## Get afv for design points
  # Not needed in this version of the function
  # afvColumn <- glmssn$args$addfunccol
  # design.afv <- ssn@obspoints@SSNPoints[[1]]@point.data[ind, afvColumn]
  # 
  ## Get distance matrices, etc.
  
  #dist.junc.obs <- SSN:::getStreamDistMatInt(ssn, design.points, "obs")
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  
  # print(mat)
  
  # This code is OBSOLETE
  # distance.matrix <- important.matrices$Distance
  # connectivity.matrix <- important.matrices$Connectivity
  # weights.vector <- rep(1, nrow(connectivity.matrix))
  # conn.weight.matrix <- connectivity.matrix 
  # 
  ## Get simulated covariance parameters
  
  cvp.cols <- length(glmssn$estimates$theta)
  cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  
  # OBSOLETE CODE
  # X <- model.matrix(object = formula, data = design)
  # Xt <- t(X)
  
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
    # print(theta.i)
    V <- SSN:::makeCovMat(
      theta.i, 
      mat$d, 
      mat$a, 
      mat$b, 
      mat$w, 
      NULL, 
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
    
    ##	if(class(qrV) != "qr") {browser()}
    covbi <- t(X) %*% solve(V) %*% X
    covb <- solve(covbi)
    D[i] <- -log(det(covb))
    
  }
  
  return(mean(D, na.rm = TRUE))
  
}
