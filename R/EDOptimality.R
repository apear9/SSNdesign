#' A utility function for optimal designs for the empirical estimation of fixed effects parameters.
#' 
#'@description
#'
#'\code{EDOoptimality} is a utility function that can be used with \code{\link{findOptimalDesign}}. It is a utility function that minimises the determinant of the variance-covariance matrix of the fixed effects over a set of designs. 
#' 
#'@usage
#'
#'\code{EDOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
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
#'This utility function is similar to \code{\link{DOptimality}}, except that it does not assume the covariance parameters are known. These are instead estimated from the data associated with a set of design points, while still incorporating prior information about the covariance parameters.  
#'
#' @export
EDOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Cut down SSN to contain only the design and prediction points
  
  ind <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  ssn@obspoints@SSNPoints[[1]]@network.point.coords <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind, ]
  ssn@obspoints@SSNPoints[[1]]@point.coords <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind, ]
  ssn@obspoints@SSNPoints[[1]]@point.data <- ssn@obspoints@SSNPoints[[1]]@point.data[ind, ]
  current.netID <- as.character(ssn@obspoints@SSNPoints[[1]]@point.data$netID)
  ssn@obspoints@SSNPoints[[1]]@point.data$netID <- as.factor(as.character(ssn@obspoints@SSNPoints[[1]]@point.data$netID))
  ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID <- as.factor(as.character(ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID))
  
  # Cut down matrices involving the observations
  
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  
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
      ssn.object = ssn,
      ObsSimDF = ssn@obspoints@SSNPoints[[1]]@point.data,
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
      matrices.obs = mat
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
      w = mat$w
    )
    
    # Obtain determinant of estimated covariance matrix on the fixed effects
    
    ED[i] <- -log(det(mdl.tmp$estimates$covb))
    
  }
  
  ED <- mean(ED)
  
  return(ED)
  
}