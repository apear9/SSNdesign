#'@inherit DOptimality
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
  fep <- extra.arguments$Empirical.FEP
  
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
      CorParms = prior.parameters[i, ],
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