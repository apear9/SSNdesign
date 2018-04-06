#' A dual utility function for optimal designs for prediction with estimated covariance and fixed effects parameters.
#' 
#'@description
#'\code{EKoptimality} is a utility function that can be used with either \code{\link{findOptimalDesign}} or \code{\link{doAdaptiveDesign}}. It is a utility function that minimises the total kriging variance over a set of prediction points given a set of design points.
#' 
#'@usage
#'
#'\code{EKOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
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
#'This function implements the K-optimal design for prediction discussed in Som et al. (2014). This utility assumes the covariance parameters are unknown to the user and must be empirically estimated. The method has been slightly modified to accommodate the Monte Carlo draws which are used to approximate the utility given prior information on the covariance parameters. See \code{\link{KOptimality}} for a modified version of this utility where the covariance parameters are known and need not be estimated from the data. 
#'  
#'@export
EKOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Check that prediction sites actually exist
  
  if(length(ssn@predpoints@SSNPoints) < 1){
    
    stop("This SSN is missing prediction points. The utility cannot be evaluated.")
    
  }
  
  # Cut down SSN to contain only the design and prediction points
  
  ind <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  ssn@obspoints@SSNPoints[[1]]@network.point.coords <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind, ]
  ssn@obspoints@SSNPoints[[1]]@point.coords <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind, ]
  ssn@obspoints@SSNPoints[[1]]@point.data <- ssn@obspoints@SSNPoints[[1]]@point.data[ind, ]
  
  current.netID <- as.character(ssn@obspoints@SSNPoints[[1]]@point.data$netID)
  ind <- ssn@predpoints@SSNPoints[[1]]@point.data$netID %in% current.netID
  ssn@predpoints@SSNPoints[[1]]@network.point.coords <- ssn@predpoints@SSNPoints[[1]]@network.point.coords[ind, ]
  ssn@predpoints@SSNPoints[[1]]@point.coords <- ssn@predpoints@SSNPoints[[1]]@point.coords[ind, ]
  ssn@predpoints@SSNPoints[[1]]@point.data <- ssn@predpoints@SSNPoints[[1]]@point.data[ind, ]
  
  # Ensure that the only level of the factor netID is the current netID
  
  ssn@obspoints@SSNPoints[[1]]@point.data$netID <- as.factor(as.character(ssn@obspoints@SSNPoints[[1]]@point.data$netID))
  ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID <- as.factor(as.character(ssn@obspoints@SSNPoints[[1]]@network.point.coords$NetworkID))
  ssn@predpoints@SSNPoints[[1]]@point.data$netID <- as.factor(as.character(ssn@predpoints@SSNPoints[[1]]@point.data$netID))
  ssn@predpoints@SSNPoints[[1]]@network.point.coords$NetworkID <- as.factor(as.character(ssn@predpoints@SSNPoints[[1]]@network.point.coords$NetworkID))
  
  # Cut down matrices involving the observations
  
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  extra.arguments$Matrices.obs <- mat
  
  # Do the same for the pxo matrix
  mat <- extra.arguments$Matrices.pxo
  mat$d <-  mat$d[ind.mat, ]
  mat$a <-  mat$a[ind.mat, ]
  mat$b <-  mat$b[ind.mat, ] 
  mat$w <-  mat$w[ind.mat, ]
  extra.arguments$Matrices.pxo <- mat
  
  # Simulate parameters as required
  
  cvp.cols <- length(glmssn$estimates$theta)
  cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  
  for(i in 1:cvp.cols){
    
    cvp[, i] <- prior.parameters[[i]](n.draws) # The covariance parameters
    
  }
  
  fep <- MASS::mvrnorm(n.draws, glmssn$estimates$betahat, glmssn$estimates$covb) # The fixed effects
  fep <- unname(fep)
  
  # Estimate utility
  
  EK <- vector("numeric", n.draws)
  
  mod.formula <- as.character(glmssn$args$formula)
  nterms <- length(mod.formula)
  mod.formula <- as.formula(mod.formula[c(1, 3:nterms)])
  
  for(i in 1:n.draws){
    
    # Simulate data from simulated FE and CovParms values
    
    print("here1")
    
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
      matrices.obs = extra.arguments$Matrices.obs,
      matrices.preds = extra.arguments$Matrices.prd,
      matrices.predsxobs = extra.arguments$Matrices.pxo
    )$ssn.object
    
    # Fit model to simulated data
    
    print("here2")
    
    mdl.tmp <- glmssn_minimal(
      formula = glmssn$args$formula,
      ssn.object = ssn.i,
      family = glmssn$args$family,
      CorModels = glmssn$args$CorModels,
      use.nugget = glmssn$args$use.nugget,
      use.anisotropy = glmssn$args$use.anisotropy,
      addfunccol = glmssn$args$addfunccol,
      useTailDownWeight = glmssn$args$useTailDownWeight,
      d = extra.arguments$Matrices.obs$d,
      a = extra.arguments$Matrices.obs$a,
      b = extra.arguments$Matrices.obs$b,
      c = extra.arguments$Matrices.obs$c,
      w = extra.arguments$Matrices.obs$w
    )
    
    # Obtain kriging variances at prediction sites from this model
    
    print("here3")
    
    n.preds.uncertainty <- ncol(glmssn$ssn.object@predpoints@SSNPoints[[1]]@point.data) + 2 # based on behaviour of predict.glmssn
    EK.i <- predict.glmssn(mdl.tmp, "preds")$ssn.object@predpoints@SSNPoints[[1]]@point.data[,n.preds.uncertainty]^2
    
    # Sum, invert
    
    EK[i] <- 1/sum(EK.i)
    
  }
  
  EK <- mean(EK)
  
  return(EK)
  
}