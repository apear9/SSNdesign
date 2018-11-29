#' A utility function for covariance parameter estimation (sequential)
#' 
#'@description
#'   
#'\code{sequentialCPOptimality} is a utility function that can be used with \code{\link{optimiseSSNDesign}}. It is a utility function that minimises the determinant of the variance-covariance matrix of the inverse Fisher information matrix over a set of design points. 
#' 
#'@usage
#'
#'\code{sequentialCPOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
#' 
#'@param ssn An object of class SpatialStreamNetwork
#'@param glmssn An model object of class glmssn.
#'@param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#'@param prior.parameters A list of random functions that are parameterised in terms of n.draws.
#'@param n.draws A numeric scalar for the number of Monte Carlo draws to use when approximating the utility. 
#'@param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{optimiseSSNDesign}}.
#'@return A numeric scalar.
#' 
#'@details
#'
#'Note that this function operates differently to \code{CPOptimality}. The function \code{CPOptimality} assumes there are no sites which have already been incorporated into a design. It computes the expected Fisher information matrix (Som et al., 2014) for each set of simulated covariance parameters. In this sequential form, the observed information matrix is extracted for the sites which are fixed in the design. Then the sites which are not fixed are used to obtain the expected Fisher information. The observed and expected matrices are added, before being reduced to a scalar value to serve as the utility.
#' 
#'@export
sequentialCPOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Extract out the Hessian (observed information) from the glmssn object
  old.information <- glmssn$optimOutput$hessian
  if(is.null(glmssn$optimOutput$hessian)){
    stop("No hessian in optim output in the glmssn object.")
  }
  
  # Get design matrix
  ind <- row.names(extra.arguments$obs.X) %in% design.points
  X <- extra.arguments$obs.X[ind, ]
  Xt <- t(X)
  ind <- row.names(extra.arguments$obs.C) %in% design.points
  cds <- extra.arguments$obs.C[ind, ]
  
  ## Get distance matrices, etc.
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  n.zero <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
  
  ## Simulate covariance parameters
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
  FIM <- vector("numeric", n.draws)
  for(i in 1:n.draws){
    
    # Get covariance matrix on the data
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
    Vi <- solve(V)
    
    ##	precompute the matrix P
    P <- Vi - Vi %*% X %*% solve(Xt %*% Vi %*% X) %*% Xt %*% Vi
    
    ## estimate partial derivatives
    
    np <- length(theta.i)
    ep <- vector("list", np)
    em <- matrix(0, nrow = np, ncol = np)
    for(j in 1:np){
      theta.j <- theta.i
      theta.j[j] <- theta.i[j] + h
      V.j <- SSN:::makeCovMat(
        theta.j, 
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
      ep[[j]] <- (V.j - V)/h
    }
    
    for(j in 1:np){
      for(k in j:np){
        I_REML <- P %*% ep[[j]] %*% P %*% ep[[k]]
        em[k, j] <- em[j, k] <- 1/2 * sum(diag(I_REML))
      }
    }
    totFIM <- em + old.information
    FIM[i] <- log(det(totFIM))
    
  }
  
  if(any(is.infinite(FIM))){
    return(-1e9)
  } else {
    return(mean(FIM, na.rm = TRUE))
  }
  
}
