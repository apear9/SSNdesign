#' @inherit DOptimality
#' @export
sequentialDOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Extract out the variance-covariance matrix for the fixed effects from the glmssn object
  old.covbi <- glmssn$estimates$covbi
  
  # Retrieve relevant design matirx
  ind <- row.names(extra.arguments$obs.X) %in% design.points
  X <- extra.arguments$obs.X[ind,]
  Xt <- t(X)
  
  # Retrieve coordinates of sampling points
  ind <- row.names(extra.arguments$obs.C) %in% design.points
  cds <- extra.arguments$obs.C[ind,]
  
  # Cut down distance matrices, etc.
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  n.zero <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
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
  
  ## Perform MC simulations
  D <- vector("numeric", n.draws)
  
  for(i in 1:n.draws){
    
    theta.i <- prior.parameters[i, ]
    V <- SSN:::makeCovMat(
      theta.i, 
      mat$d, 
      mat$a, 
      mat$b, 
      mat$w, 
      n.zero, 
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
  #  covb <- solve(covbi)
    D[i] <- log(det(covbi + old.covbi))
    
  }
  
  return(mean(D, na.rm = TRUE))
  
}
