#' A utility function for prediction with known covariance parameters
#'
#'@description 
#'   
#' The function \code{KOptimality} implements the K-optimal utility function from Falk et al. (2014).
#' 
#'@usage
#'
#'\code{KOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)}
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
#' \deqn{U(d,\theta) = \left(\sum_{i = 1}^n Var(\hat{y}_i)\right)^{-1}}{U(d,\theta)= (\sum_{i = 1}^n Var(\hat{y}_i))^{-1}}
#' That is, it is the inverse sum of the estimated kriging variances for a set of prediction sites with indices \eqn{i = 1, 2, ..., n}. 
#' 
#' One final note: do not worry about passing arguments to this function. All arguments are dealt with internally by \code{\link{optimiseSSNDesign}}.
#'
#' @examples 
#'  
#'\dontrun{
#' # Create stream network
#' s <- createSSN(100, systematicDesign(0.25), systematicDesign(0.25), paste(tempdir(), "s.ssn", sep = "/"), TRUE)
#' createDistMat(s, "preds", TRUE, TRUE)
#' 
#' # Simulate data on network
#' s <- SimulateOnSSN(s, getSSNdata.frame(s), getSSNdata.frame(s, "preds"), "preds", formula = ~1, coefficients = 1, CorParms = c(1,2,1,2,1,2,0.1),addfunccol = "addfunccol")$ssn.object
#' 
#' # Fit a model to the simulated data
#' m <- glmssn(Sim_Values ~ 1, s, addfunccol = "addfunccol")
#' 
#' # Define the priors on the covariance parameters
#' p <- constructLogNormalCovPriors(m)
#' 
#' # Find the optimal design using D-optimality as the utility function
#' r <- optimiseSSNDesign(s, paste(tempdir(), "r.ssn", sep = "/"), m, 25, utility.function = KOptimality, prior.parameters = p)
#' 
#' # Plot result to check
#' plot(r$ssn.new, "Sim_Values")
#' }
#'      
#' @export
KOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  ## Get info for the obs sites  
  ind <- row.names(extra.arguments$obs.X) %in% design.points
  X <- extra.arguments$obs.X[ind, ]
  Xt <- t(X)
  
  ## Get coordinates
  ind.cds <- row.names(extra.arguments$obs.C) %in% design.points
  cds.obs <- extra.arguments$obs.C[ind.cds, ]
  
  ## Cut down matrices involving the observations
  
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  extra.arguments$Matrices.Obs <- mat
  net.zero.obs <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
  
  # Do the same for the pxo matrix
  mat <- extra.arguments$Matrices.pxo
  mat$d <-  mat$d[ind.mat, ]
  mat$a <-  mat$a[ind.mat, ]
  mat$b <-  mat$b[ind.mat, ] 
  mat$w <-  mat$w[ind.mat, ]
  extra.arguments$Matrices.pxo <- mat
  net.zero.pxo <- extra.arguments$net.zero.pxo[ind.mat, ]
  # Get info for the pred sites
  indp <- row.names(extra.arguments$prd.X) %in% ssn@predpoints@SSNPoints[[1]]@point.data$pid
  X0 <- t(extra.arguments$prd.X[indp, ]) # have this extracted and put in extra.arguments
  cds.prd <- extra.arguments$prd.C[indp, ]

  # Matrices
  
  ## Get other model parameters
  td <- glmssn$args$useTailDownWeight
  cm <- glmssn$args$CorModels
  un <- glmssn$args$use.nugget
  ua <- glmssn$args$use.anisotropy
  re <- glmssn$sampInfo$REs
  
  # ## Simulate covparms from priors
  # cvp.cols <- length(glmssn$estimates$theta)
  # cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  # for(i in 1:cvp.cols){
  #   cvp[, i] <- prior.parameters[[i]](n.draws)
  # }
  
  ## Loop here to find utilities across all simulations
  
  # initialise empty vector to store results
  
  K_all <- vector("numeric", n.draws)
  
  for(i in 1:n.draws){
    
    # get simulated covariance parameters
    
    theta.i <- prior.parameters[i,]
    
    # get W and Wi
    
    W <- SSN:::makeCovMat(
      theta.i, 
      extra.arguments$Matrices.Obs$d, 
      extra.arguments$Matrices.Obs$a, 
      extra.arguments$Matrices.Obs$b, 
      extra.arguments$Matrices.Obs$w, 
      net.zero.obs, 
      cds.obs[, "x"], 
      cds.obs[, "y"], 
      cds.obs[, "x"], 
      cds.obs[, "y"], 
      td,
      cm, 
      un,
      ua,
      re
    )
    Wi <- solve(W)
    
    # get V
    
    V <- SSN:::makeCovMat(
      theta.i, 
      extra.arguments$Matrices.prd$d, 
      extra.arguments$Matrices.prd$a, 
      extra.arguments$Matrices.prd$b, 
      extra.arguments$Matrices.prd$w, 
      extra.arguments$net.zero.prd, 
      cds.prd[, "x"], 
      cds.prd[, "y"], 
      cds.prd[, "x"], 
      cds.prd[, "y"], 
      td,
      cm, 
      un,
      ua,
      re
    )
    
    # get C and tC
    
    C <- SSN:::makeCovMat(
      theta.i, 
      extra.arguments$Matrices.pxo$d, 
      extra.arguments$Matrices.pxo$a, 
      extra.arguments$Matrices.pxo$b, 
      extra.arguments$Matrices.pxo$w, 
      net.zero.pxo, 
      cds.obs[, "x"], 
      cds.obs[, "y"], 
      cds.prd[, "x"], 
      cds.prd[, "y"], 
      td,
      cm, 
      FALSE,
      ua,
      re
    )
    Ct <- t(C)
    
    # get M and Mt
    
    M <- (X0 - Xt %*% Wi %*% C)
    Mt <- t(M)
    
    # pred utility
    
    K <- tryCatch(V - Ct %*% Wi %*% C + Mt %*% solve(Xt %*% Wi %*% X) %*% M, error = function(e) return(matrix(c(-1e9,0,0,0), nrow = 2)))
    
    # get the trace of the matrix at the observed sites
    
    K_all[i] <- 1/sum(diag(K))
    
  }
  
  # # spit out result
  Ud <- mean(K_all)
  if(Ud <= 0){ # Change to any(K_all)
    return(-1e9)
  } else {
    return(Ud)
  } 
  
}
