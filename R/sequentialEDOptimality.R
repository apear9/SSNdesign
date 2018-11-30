#' A utility function for empirical parameter estimation (sequential)
#'
#'@description
#'   
#'\code{sequentialEDOptimality} is a utility function that can be used with \code{\link{optimiseSSNDesign}}. It is a utility function that minimises the determinant of the variance-covariance matrix of the fixed effects for models that are fit over a set of design points. 
#' 
#'@usage
#'
#'\code{sequentialEDOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
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
#'Note that this function operates differently to \code{EDOptimality}. The functions \code{DOptimality} and \code{EDOptimality} assume there are no sites which have already been incorporated into a design. They compute the variance-covariance matrix on the fixed effects (Som et al., 2014) for each set of simulated or esimated covariance parameters, respectively. In this sequential form, the observed variance-covariance matrix is extracted for the sites which are fixed in the design. Then the sites which are not fixed are used to estimate the expected variance-covariance matrix. The observed and expected matrices are added, before being reduced to a scalar value to serve as the utility.
#' 
#'@export
sequentialEDOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Extract out the variance-covariance matrix for the fixed effects from the glmssn object
  old.covbi <- glmssn$estimates$covbi
  
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
  n.zero <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
  
  # Simulate parameters as required
  cvp.cols <- length(glmssn$estimates$theta)
  cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  for(i in 1:cvp.cols){
    cvp[, i] <- prior.parameters[[i]](n.draws) # The covariance parameters
  }
  fep <- MASS::mvrnorm(n.draws, glmssn$estimates$betahat, glmssn$estimates$covb) # The fixed effects
  fep <- unname(fep)
  
  # Extract model formula
  mod.formula <- as.character(glmssn$args$formula)
  nterms <- length(mod.formula)
  mod.formula <- as.formula(mod.formula[c(1, 3:nterms)])
  
  # Estimate utility
  ED <- vector("numeric", n.draws)
  
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
    ED[i] <- log(det(mdl.tmp$estimates$covbi + old.covbi))
    
  }
  
  ED <- mean(ED)
  
  return(ED)
  
}