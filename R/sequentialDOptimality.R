#' A utility function for parameter estimation with known covariance parameters (sequential)
#' 
#' \code{sequentialDOptimality()}
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param glmssn An object of class glm.ssn
#' @param design.points A numeric vector of pids for sites in a proposed design
#' @param prior.parameters A list of functions parameterised in terms of n.draws that simulate from probability distributions
#' @param extra.arguments A list of extra arguments that the utility may rely on
#' @return A numeric representing the sequential CP-optimal utility of a design estimated by Monte-Carlo integration
#' 
#' @export
sequentialDOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Extract out the variance-covariance matrix for the fixed effects from the glmssn object
  obs.covb <- glmssn$estimates$covb
  
  # Exclude the fixed points from the set of design points
  ind <- !(design.points %in% extra.arguments$fixed)
  design.points <- design.points[ind]
  
  # Retrieve relevant design matirx
  ind <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  X <- glmssn$sampinfo$X[ind, ]
  Xt <- t(X)
  
  # Retrieve coordinates of sampling points
  cds <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind, ]
  colnames(cds) <- c("x", "y")
  
  # Cut down distance matrices, etc.
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  
  ## Get other model parameters
  td <- glmssn$args$useTailDownWeight
  cm <- glmssn$args$CorModels
  un <- glmssn$args$use.nugget
  ua <- glmssn$args$use.anisotropy
  re <- glmssn$sampInfo$REs
  
  ## Simulate covparms from priors
  cvp.cols <- length(glmssn$estimates$theta)
  cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  for(i in 1:cvp.cols){
    cvp[, i] <- prior.parameters[[i]](n.draws)
  }
  
  ## Perform MC simulations
  D <- vector("numeric", n.draws)
  
  for(i in 1:n.draws){
    
    theta.i <- cvp[i, ]
    V <- SSN:::makeCovMat(
      theta.i, 
      mat$d, 
      mat$a, 
      mat$b, 
      mat$w, 
      NULL, 
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
    
    covbi <- t(X) %*% solve(V) %*% X
    covb <- solve(covbi)
    D[i] <- -log(det(covb + obs.covb))
    
  }
  
  return(mean(D, na.rm = TRUE))
  
}
