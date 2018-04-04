KOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Get info for the obs sites  
  ind <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  X <- glmssn$sampinfo$X[ind, ]
  Xt <- t(X)
  cds.obs <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind, ]
  colnames(cds.obs) <- c("x", "y")
  
  # Cut down matrices involving the observations
  
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  extra.arguments$Matrices.Obs <- mat
  
  # Do the same for the pxo matrix
  
  mat <- extra.arguments$Matrices.pxo
  mat$d <-  mat$d[ind.mat, ]
  mat$a <-  mat$a[ind.mat, ]
  mat$b <-  mat$b[ind.mat, ] 
  mat$w <-  mat$w[ind.mat, ]
  extra.arguments$Matrices.pxo <- mat
  
  # Get info for the pred sites
  this.net <- ssn@obspoints@SSNPoints[[1]]@point.data$netID[ind][1]
  ind <- ssn@predpoints@SSNPoints[[1]]@point.data$netID == this.net
  X0 <- t(model.matrix(glmssn$args$formula, ssn@predpoints@SSNPoints[[1]]@point.data[ind,]))
  cds.prd <- ssn@predpoints@SSNPoints[[1]]@point.coords[ind, ]
  colnames(cds.prd) <- c("x", "y")
  # Matrices
  
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
  
  ## Loop here to find utilities across all simulations
  
  # initialise empty vector to store results
  
  K_all <- vector("numeric", n.draws)

  for(i in 1:n.draws){
    
    # get simulated covariance parameters
    
    theta.i <- cvp[i, ]
    
    # get W and Wi
    
    W <- SSN:::makeCovMat(
      theta.i, 
      extra.arguments$Matrices.Obs$d, 
      extra.arguments$Matrices.Obs$a, 
      extra.arguments$Matrices.Obs$b, 
      extra.arguments$Matrices.Obs$w, 
      NULL, 
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
      NULL, 
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
      NULL, 
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
    
    K <- V - Ct %*% Wi %*% C + Mt %*% solve(Xt %*% Wi %*% X) %*% M
    
    # get the trace of the matrix at the observed sites
    
    K_all[i] <- 1/sum(diag(K))
    
  }
  
  # spit out result
  
  return(mean(K_all))
  
}