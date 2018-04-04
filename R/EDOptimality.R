#' A utility function for optimal designs for the empirical estimation of fixed effects parameters.
#' 
#' \code{EDOptimality()}
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param glmssn An object of class glm.ssn
#' @param design.points A numeric vector of pids for sites in a proposed design
#' @param prior.parameters A list of functions parameterised in terms of n.draws that simulate from probability distributions
#' @param extra.arguments A list of extra arguments that the utility may rely on
#' @return A numeric representing the D-optimal utility of a design estimated by Monte-Carlo integration
#' 
#' @export
EDOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
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