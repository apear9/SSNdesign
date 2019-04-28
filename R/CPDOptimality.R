#' A dual-purpose utility function for covariance and fixed effect parameter estimation
#' 
#'@description
#'   
#' The function \code{CPDOptimality} combines the D- and CPP-optimal utility functions from Falk et al. (2014).
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param glmssn A fitted model object of class glmssn.
#' @param design.points A vector of pids corresponding to a set of observed sites in the obspoints slot of the SpatialStreamNetwork object.
#' @param prior.parameters A list of random functions that are parameterised in terms of n.draws.
#' @param n.draws A numeric scalar for the number of Monte Carlo draws to use when approximating the utility. 
#' @param extra.arguments A list of extra parameters that control the behaviour of the utility function. The distance matrices required to compute covariance matrices are also stored in this list. Note that these are generated inside \code{\link{optimiseSSNDesign}}.
#' @return A single number representing the expected utility for the design specified by \code{design.points}.
#' 
#' @details A value of \code{-1e9} or \code{-2e9} indicates a failure of this utility function. Usually this means either the design matrix or the covariance matrix is rank-deficient. 
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
#' r <- optimiseSSNDesign(s, paste(tempdir(), "r.ssn", sep = "/"), m, 25, utility.function = CPDOptimality, prior.parameters = p)
#' 
#' # Plot result to check
#' plot(r$ssn.new, "Sim_Values")
#' }
#' 
#' @export
CPDOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  CPOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments) + DOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)
  
}