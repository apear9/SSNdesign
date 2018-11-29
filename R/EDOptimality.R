#' A utility function for optimal designs for the empirical estimation of fixed effects parameters.
#' 
#'@description
#'
#' The function \code{EDOoptimality} implements the ED-optimal utility function from Falk et al. (2014).
#' 
#'@usage
#'
#'\code{EDOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
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
#' \deqn{U(d, \theta, y) = \log\det(\hat{Var}(\beta)^{-1})}{U(d, \theta, y) = log(det(\hat{Var}(\beta)^{-1}))}
#' 
#' where \eqn{\hat{Var}(\beta)} is the estimated variance-covariance matrix of the fixed effects from a fitted \code{glmssn} object.
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
#' r <- optimiseSSNDesign(s, paste(tempdir(), "r.ssn", sep = "/"), m, 25, utility.function = EDOptimality, prior.parameters = p)
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
EDOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Cut down SSN to contain only the design and prediction points
  ssn2 <- ssn
  ind.x <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  ind.c <- row.names(ssn@obspoints@SSNPoints[[1]]@point.coords) %in% design.points
  ind.n <- row.names(ssn@obspoints@SSNPoints[[1]]@network.point.coords) %in% design.points
  ssn2@obspoints@SSNPoints[[1]]@point.data <- ssn2@obspoints@SSNPoints[[1]]@point.data[ind.x, ]
  ssn2@obspoints@SSNPoints[[1]]@point.coords <- ssn2@obspoints@SSNPoints[[1]]@point.coords[ind.c, ]
  ssn2@obspoints@SSNPoints[[1]]@network.point.coords <- ssn2@obspoints@SSNPoints[[1]]@network.point.coords[ind.n, ]
  
  # Cut down matrices involving the observations
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  n.zero <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
  
  # Simulate parameters as required
  cvp.cols <- length(glmssn$estimates$theta)
  cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  for(i in 1:cvp.cols){
    cvp[, i] <- prior.parameters[[i]](n.draws) # The covariance parameters
  }
  fep <- MASS::mvrnorm(n.draws, glmssn$estimates$betahat, glmssn$estimates$covb) # The fixed effects
  fep <- unname(fep)
  
  # Estimate utility
  ED <- vector("numeric", n.draws)
  mod.formula <- as.character(glmssn$args$formula)
  nterms <- length(mod.formula)
  mod.formula <- as.formula(mod.formula[c(1, 3:nterms)])
  for(i in 1:n.draws){
    # Simulate data from simulated FE and CovParms values
    ssn.i <- SimulateOnSSN_minimal( 
      ssn.object = ssn2,
      ObsSimDF = ssn2@obspoints@SSNPoints[[1]]@point.data,
      PredSimDF = NULL,
      PredID = NULL,
      formula = mod.formula,
      coefficients = fep[i, ],
      CorModels = glmssn$args$CorModels,
      use.nugget = glmssn$args$use.nugget,
      CorParms = cvp[i, ],
      use.anisotropy = glmssn$args$use.anisotropy,
      addfunccol = glmssn$args$addfunccol,
      useTailDownWeight = glmssn$args$useTailDownWeight,
      family = glmssn$args$family,
      matrices.obs = mat,
      net.zero.obs = n.zero
    )$ssn.object
    # Fit model to simulated data
    mdl.tmp <- glmssn_minimal( 
      formula = glmssn$args$formula,
      ssn.object = ssn.i,
      family = glmssn$args$family,
      CorModels = glmssn$args$CorModels,
      use.nugget = glmssn$args$use.nugget,
      use.anisotropy = glmssn$args$use.anisotropy,
      addfunccol = glmssn$args$addfunccol,
      useTailDownWeight = glmssn$args$useTailDownWeight,
      d = mat$d,
      a = mat$a,
      b = mat$b,
      c = mat$c,
      w = mat$w,
      n = n.zero
    )
    # Obtain determinant of estimated covariance matrix on the fixed effects
    ED[i] <- log(det(mdl.tmp$estimates$covbi))
  }
  
  # Catch errors
  if(any(is.infinite(ED))){
    ED <- -1e9
  } else {
    ED <- mean(ED)
  }
  return(ED)
  
}
