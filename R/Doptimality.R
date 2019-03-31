#' A utility function for fixed effects parameter estimation with known covariance parameters
#' 
#' @description
#'
#' The function \code{DOoptimality} implements the D-optimal utility function from Falk et al. (2014). 
#' 
#' This function can be used with \code{\link{optimiseSSNDesign}}. 
#' 
#' @usage
#'
#' \code{DOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param glmssn A fitted model object of class glmssn.
#' @param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#' @param prior.parameters A list of random functions that are parameterised in terms of n.draws.
#' @param n.draws A numeric scalar for the number of Monte Carlo draws to use when approximating the utility. 
#' @param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{optimiseSSNDesign}}.
#' @return A single number representing the expected utility for the design specified by \code{design.points}.
#' 
#' @details
#'
#' The utility function is
#' 
#' \deqn{U(d, \theta) = \log\det(X'\Sigma(\theta)^-1X)}{U(d, \theta) = log(det(X'\Sigma(\theta)^-1X))}
#' 
#' where \eqn{X} is the design matrix for the model specified by \code{glmssn}, \eqn{X'} is its transpose, and \eqn{\Sigma(\theta)} is the covariance matrix for the data, which is obtained by evaluating a covariance function given a distance matrix and simulated values for the parameters \eqn{\theta}. The values of \eqn{\theta} are simulated from the priors in \code{prior.parameters}. The result will be averaged over \code{n.draws} simulations of the parameters \eqn{\theta}.
#' 
#' One final note: do not worry about passing arguments to this function. All arguments are handled internally by \code{\link{optimiseSSNDesign}}. 
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#' # Create stream network
#' s <- createSSN(100, systematicDesign(0.25), path = paste(tempdir(), "s.ssn", sep = "/"), importToR = TRUE)
#' createDistMat(s)
#' 
#' # Simulate data on network
#' s <- SimulateOnSSN(s, getSSNdata.frame(s), formula = ~1, coefficients = 1, CorParms = c(1,2,1,2,1,2,0.1),addfunccol = "addfunccol")$ssn.object
#' 
#' # Fit a model to the simulated data
#' m <- glmssn(Sim_Values ~ 1, s, addfunccol = "addfunccol")
#' 
#' # Define the priors on the covariance parameters
#' p <- constructLogNormalCovPriors(m)
#' 
#' # Find the optimal design using D-optimality as the utility function
#' r <- optimiseSSNDesign(s, paste(tempdir(), "r.ssn", sep = "/"), m, 25, utility.function = DOptimality, prior.parameters = p)
#' 
#' # Plot result to check
#' plot(r$ssn.new, "Sim_Values")
#' 
#' }
#' 
#' @references 
#' 
#' Falk, M., McGree, J.M., and Pettit, A.N. (2014). Sampling designs on stream networks using the pseudo-Bayesian approach. \emph{Environmental and Ecological Statistics}, \emph{21}(4), 751-773.
#' 
#' @export
DOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  #t1 <- Sys.time()
  
  ## Get data for design points. Note that design matrix is stored as obs.X in extra.arguments
  ind <- row.names(extra.arguments$obs.X) %in% design.points 
  X <- extra.arguments$obs.X[ind, ]
  Xt <- t(X)
  
  ## Get coordinates of points in case the covariance function includes a Euclidean dist-based component
  #  Note that the coordinates (with the same ordering as the design matrix) is stored as obs.C in extra.arguments
  ind.cds <- row.names(extra.arguments$obs.C) %in% design.points
  cds <- extra.arguments$obs.C[ind.cds, ] # names are already x and y
  
  ## Get distance matrices, etc.
  
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
  
  ## Perform MC simulations
  
  D <- vector("numeric", n.draws)
  
  for(i in 1:n.draws){
    
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
    covbi <- Xt %*% solve(V) %*% X
    
    D[i] <- log(det(covbi)) 
    # Note that this is the same as -log(det(solve(covbi))) since the inverse of a determinant is the determinant of the inverse
    # Losing the inversion saves time, especially for large numbers of Monte Carlo draws.
    
  }
  
  if(any(is.infinite(D))){
    return(-1e9)
  } else {
    return(mean(D, na.rm = TRUE))
  }
  
}
