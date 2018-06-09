#' A utility function for prediction with known covariance parameters
#'
#'@description 
#'   
#'\code{Koptimality} is a utility function that can be used with either \code{\link{findOptimalDesign}} or \code{\link{doAdaptiveDesign}}. It is a utility function that minimises the total kriging variance over a set of prediction points given a set of design points.
#' 
#'@usage
#'
#'\code{KOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
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
#'This function implements the K-optimal designs for prediction discussed in Som et al. (2014). This utility assumes the covariance parameters are known to the user and do not have to be empirically estimated. The method has been slightly modified to accommodate the Monte Carlo draws which are used to approximate the utility given prior information on the covariance parameters. See \code{\link{EKOptimality}} for a modified version of this utility where the covariance parameters are unknown and must be estimated from the data. 
#'  
#'@export
KOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  #t1 <- Sys.time()
  # Get info for the obs sites  
  ind <- row.names(glmssn$sampinfo$X) %in% design.points
  X <- glmssn$sampinfo$X[ind, ]
  Xt <- t(X)
  # Coordinates
  ind.cds <- row.names(ssn@obspoints@SSNPoints[[1]]@point.coords) %in% design.points
  cds.obs <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind.cds, ]
  colnames(cds.obs) <- c("x", "y")
  
  # Cut down matrices involving the observations
  
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  extra.arguments$Matrices.Obs <- mat
  net.zero.obs <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
  
  # Do the same for the pxo matrix
  mat <- extra.arguments$Matrices.pxo
  mat$d <-  mat$d[ind.mat, ]
  mat$a <-  mat$a[ind.mat, ]
  mat$b <-  mat$b[ind.mat, ] 
  mat$w <-  mat$w[ind.mat, ]
  extra.arguments$Matrices.pxo <- mat
  net.zero.pxo <- extra.arguments$net.zero.pxo[ind.mat, ]
    
  # Get info for the pred sites
  # this.net <- ssn@obspoints@SSNPoints[[1]]@point.data$netID[ind][1]
  # ind <- ssn@predpoints@SSNPoints[[1]]@point.data$netID == this.net
  indp <- ssn@predpoints@SSNPoints[[1]]@point.data$pid %in% row.names(extra.arguments$Matrices.prd$d)
  X0 <- t(model.matrix(glmssn$args$formula, ssn@predpoints@SSNPoints[[1]]@point.data[indp, ]))
  cds.prd <- ssn@predpoints@SSNPoints[[1]]@point.coords[indp, ]
  colnames(cds.prd) <- c("x", "y")
  # Matrices
  
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
  
  ## Loop here to find utilities across all simulations
  
  # initialise empty vector to store results
  
  K_all <- vector("numeric", n.draws)

  for(i in 1:n.draws){
    
    # get simulated covariance parameters
    
    theta.i <- cvp[i, ]
    
    # get W and Wi
    
    W <- SSN:::makeCovMat(
      theta.i, 
      extra.arguments$Matrices.Obs$d, 
      extra.arguments$Matrices.Obs$a, 
      extra.arguments$Matrices.Obs$b, 
      extra.arguments$Matrices.Obs$w, 
      net.zero.obs, 
      cds.obs[, "x"], 
      cds.obs[, "y"], 
      cds.obs[, "x"], 
      cds.obs[, "y"], 
      td,
      cm, 
      un,
      ua,
      re
    )
    Wi <- solve(W)
    
    # get V
    
    V <- SSN:::makeCovMat(
      theta.i, 
      extra.arguments$Matrices.prd$d, 
      extra.arguments$Matrices.prd$a, 
      extra.arguments$Matrices.prd$b, 
      extra.arguments$Matrices.prd$w, 
      extra.arguments$net.zero.prd, 
      cds.prd[, "x"], 
      cds.prd[, "y"], 
      cds.prd[, "x"], 
      cds.prd[, "y"], 
      td,
      cm, 
      un,
      ua,
      re
    )
    
    # get C and tC
    
    C <- SSN:::makeCovMat(
      theta.i, 
      extra.arguments$Matrices.pxo$d, 
      extra.arguments$Matrices.pxo$a, 
      extra.arguments$Matrices.pxo$b, 
      extra.arguments$Matrices.pxo$w, 
      net.zero.pxo, 
      cds.obs[, "x"], 
      cds.obs[, "y"], 
      cds.prd[, "x"], 
      cds.prd[, "y"], 
      td,
      cm, 
      FALSE,
      ua,
      re
    )
    Ct <- t(C)
    
    # get M and Mt
    
    M <- (X0 - Xt %*% Wi %*% C)
    Mt <- t(M)
    
    # pred utility
    
    K <- V - Ct %*% Wi %*% C + Mt %*% solve(Xt %*% Wi %*% X) %*% M
    
    # get the trace of the matrix at the observed sites
    
    K_all[i] <- 1/sum(diag(K))
    
  }
  
  # # spit out result
  #print(Sys.time() - t1)
  return(mean(K_all))
  
}
