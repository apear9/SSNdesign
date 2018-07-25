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
  # t1 <- Sys.time()
  # Cut down SSN to contain only the design and prediction points
  
  ind <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  ssn@obspoints@SSNPoints[[1]]@network.point.coords <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind, ]
  ssn@obspoints@SSNPoints[[1]]@point.coords <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind, ]
  ssn@obspoints@SSNPoints[[1]]@point.data <- ssn@obspoints@SSNPoints[[1]]@point.data[ind, ]
  
  # Get model matrix anyway
  
  ind <- row.names(glmssn$sampinfo$X) %in% design.points
  X <- glmssn$sampinfo$X[ind, ]
  Xt <- t(X)
  cds.obs <- ssn@obspoints@SSNPoints[[1]]@point.coords
  colnames(cds.obs) <- c("x", "y")
  
  # Get info for the pred sites
  # this.net <- ssn@obspoints@SSNPoints[[1]]@point.data$netID[ind][1]
  # ind <- ssn@predpoints@SSNPoints[[1]]@point.data$netID == this.net
  
  indp <- ssn@predpoints@SSNPoints[[1]]@point.data$pid %in% row.names(extra.arguments$Matrices.prd$d)
  # X0 <- t(model.matrix(glmssn$args$formula, ssn@predpoints@SSNPoints[[1]]@point.data[indp, ]))
  cds.prd <- ssn@predpoints@SSNPoints[[1]]@point.coords[indp, ]
  colnames(cds.prd) <- c("x", "y")
  
  ## Get other model parameters
  
  td <- glmssn$args$useTailDownWeight
  cm <- glmssn$args$CorModels
  un <- glmssn$args$use.nugget
  ua <- glmssn$args$use.anisotropy
  re <- glmssn$sampInfo$REs
  
  # Cut down matrices involving the observations
  
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  extra.arguments$Matrices.obs <- mat
  net.zero.obs <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
  
  # Do the same for the pxo matrix
  mat <- extra.arguments$Matrices.pxo
  mat$d <-  mat$d[ind.mat, ]
  mat$a <-  mat$a[ind.mat, ]
  mat$b <-  mat$b[ind.mat, ] 
  mat$w <-  mat$w[ind.mat, ]
  extra.arguments$Matrices.pxo <- mat
  net.zero.pxo <- extra.arguments$net.zero.pxo[ind.mat, ]
  
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
    
 #   print("here1")
    
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
      net.zero.obs = net.zero.obs,
      matrices.preds = extra.arguments$Matrices.prd,
      net.zero.preds = extra.arguments$net.zero.prd,
      matrices.predsxobs = extra.arguments$Matrices.pxo,
      net.zero.predsxobs = net.zero.pxo
    )$ssn.object
    
    # Fit model to simulated data
    
   # print("here2")
    
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
      w = extra.arguments$Matrices.obs$w,
      n = net.zero.obs
    )
    
    # Obtain kriging variances at prediction sites from this model
    
    n.preds.uncertainty <- ncol(glmssn$ssn.object@predpoints@SSNPoints[[1]]@point.data) + 2 # based on behaviour of predict.glmssn
    EK.i <- predict.glmssn(mdl.tmp, "preds")$ssn.object@predpoints@SSNPoints[[1]]@point.data[,n.preds.uncertainty]^2
    # # get estimated covariance matrix on the data
    # 
    # Wi <- mdl.tmp$estimates$Vi
    
    # get V
    
#    print("here3")
    
    # V <- SSN:::makeCovMat(
    #   as.vector(mdl.tmp$estimates$theta), 
    #   extra.arguments$Matrices.prd$d, 
    #   extra.arguments$Matrices.prd$a, 
    #   extra.arguments$Matrices.prd$b, 
    #   extra.arguments$Matrices.prd$w, 
    #   extra.arguments$net.zero.prd, 
    #   cds.prd[, "x"], 
    #   cds.prd[, "y"], 
    #   cds.prd[, "x"], 
    #   cds.prd[, "y"], 
    #   td,
    #   cm, 
    #   un,
    #   ua,
    #   re
    # )
    # 
    # # get C and tC
    # 
    # C <- SSN:::makeCovMat(
    #   as.vector(mdl.tmp$estimates$theta), 
    #   extra.arguments$Matrices.pxo$d, 
    #   extra.arguments$Matrices.pxo$a, 
    #   extra.arguments$Matrices.pxo$b, 
    #   extra.arguments$Matrices.pxo$w, 
    #   net.zero.pxo, 
    #   cds.obs[, "x"], 
    #   cds.obs[, "y"], 
    #   cds.prd[, "x"], 
    #   cds.prd[, "y"], 
    #   td,
    #   cm, 
    #   FALSE,
    #   ua,
    #   re
    # )
    # Ct <- t(C)
    # 
    # # get M and Mt
    # 
    # M <- (X0 - Xt %*% Wi %*% C)
    # Mt <- t(M)
    # 
    # # pred utility
    # 
    # EK.i <- V - Ct %*% Wi %*% C + Mt %*% solve(Xt %*% Wi %*% X) %*% M
    
    # Sum, invert
    
    EK[i] <- 1/sum(EK.i)
    
  }
  
  EK <- mean(EK)
  #print(t1 - Sys.time()) #50 out of 80 sites, with 500 draws ~~ 12 minutes
  return(EK)
  
}
