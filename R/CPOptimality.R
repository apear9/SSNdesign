#' A utility function for covariance parameter estimation
#' 
#'@description
#'   
#' The function \code{CPOptimality} implements the CP-optimal utility function from Falk et al. (2014).
#' 
#'@usage
#'
#'\code{CPOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param glmssn A fitted model object of class glmssn.
#' @param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#' @param prior.parameters A list of random functions that are parameterised in terms of n.draws.
#' @param n.draws A numeric scalar for the number of Monte Carlo draws to use when approximating the utility. 
#' @param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{optimiseSSNDesign}}.
#' @return A single number representing the expected utility for the design specified by \code{design.points}.
#' 
#'@details
#'
#' The utility function is
#' 
#' \deqn{U(d,\theta) = \log \det \left(I(\theta)\right)}{U(d, \theta) = log(det(I(\theta)))}
#' 
#' where \eqn{I(\theta)} is the expected Fisher information matrix for the covariance parameters \eqn{\theta}.
#' 
#' One final note: do not worry about passing arguments to this function. All arguments are dealt with internally by \code{\link{optimiseSSNDesign}}.
#' 
#' @examples
#' 
#' \dontrun{
#' # Create stream network
#' s <- createSSN(100, systematicDesign(0.25), path = paste(tempdir(), "s.ssn", sep = "/"), importToR = TRUE)
#' createDistMat(s)
#' 
#' # Simulate data on network
#' s <- SimulateOnSSN(s, getSSNdata.frame(s),formula = ~1, coefficients = 1, CorParms = c(1,2,1,2,1,2,0.1),addfunccol = "addfunccol")$ssn.object
#' 
#' # Fit a model to the simulated data
#' m <- glmssn(Sim_Values ~ 1, s, addfunccol = "addfunccol")
#' 
#' # Define the priors on the covariance parameters
#' p <- constructLogNormalCovPriors(m)
#' 
#' # Find the optimal design using D-optimality as the utility function
#' r <- optimiseSSNDesign(s, paste(tempdir(), "r.ssn", sep = "/"), m, 25, utility.function = CPOptimality, prior.parameters = p)
#' 
#' # Plot result to check
#' plot(r$ssn.new, "Sim_Values")
#' }
#' 
#' @export
CPOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Get design matrix
  ind <- row.names(extra.arguments$obs.X) %in% design.points
  X <- extra.arguments$obs.X[ind, ]
  Xt <- t(X)
  
  # Get coord information
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
  # cvp.cols <- length(glmssn$estimates$theta)
  # cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  
  ## Get other model parameters
  td <- glmssn$args$useTailDownWeight
  cm <- glmssn$args$CorModels
  un <- glmssn$args$use.nugget
  ua <- glmssn$args$use.anisotropy
  re <- glmssn$sampInfo$REs
  
  ## Simulate covparms from priors
  # for(i in 1:cvp.cols){
  #   cvp[, i] <- prior.parameters[[i]](n.draws)
  # }
  
  ## Perform MC simulations
  FIM <- vector("numeric", n.draws)
  
  # Before we enter the loop, make sure the step size for ffd is defined
  if(is.null(extra.arguments$h)){
    h <- 1e-5
  } else {
    h <- extra.arguments$h
  }
  
  # LOOPPPPP
  for(i in 1:n.draws){
    
    # Get covariance matrix on the data
    theta.i <- prior.parameters[i, ]
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
    
    #	precompute the matrix P
    P <- Vi - Vi %*% X %*% solve(Xt %*% Vi %*% X) %*% Xt %*% Vi
    
    # estimate partial derivatives
    np <- length(theta.i) # np = number of parameters
    ep <- vector("list", np) # empty vector of partial derivatives (hence e p)
    em <- matrix(0, nrow = np, ncol = np) # empty vector of matrices (hence e m)
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
    FIM[i] <- log(det(em)) 
    
  }
  
  if(any(is.infinite(FIM))){
    return(-1e9)
  } else{
    return(mean(FIM, na.rm = TRUE))
  }
  
}
