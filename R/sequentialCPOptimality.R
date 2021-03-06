#'@inherit DOptimality
#'@export
sequentialCPOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Get increment for numerical differentiation
  h <- extra.arguments$h 
  if(is.null(h)) h <- 1e-8 # safe default value
  
  # Extract out the Hessian (observed information) from the glmssn object
  old.information <- getCalculatedCovMatrix(glmssn)
  if(is.null(glmssn$optimOutput$hessian)){
    # Exit with error code if no hessian
    return(-1e9)
  }
  
  # Get design matrix
  ind <- row.names(extra.arguments$obs.X) %in% design.points
  X <- extra.arguments$obs.X[ind, ]
  Xt <- t(X)
  ind <- row.names(extra.arguments$obs.C) %in% design.points
  cds <- extra.arguments$obs.C[ind, ]
  
  ## Get distance matrices, etc.
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
  
  ## Perform MC simulations
  FIM <- vector("numeric", n.draws)
  for(i in 1:n.draws){
    
    # Get covariance matrix on the data
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
    Vi <- solve(V)
    
    ##	precompute the matrix P
    P <- Vi - Vi %*% X %*% solve(Xt %*% Vi %*% X) %*% Xt %*% Vi
    
    ## estimate partial derivatives
    
    np <- length(theta.i)
    ep <- vector("list", np)
    em <- matrix(0, nrow = np, ncol = np)
    for(j in 1:np){
      theta.j <- theta.i
      theta.j[j] <- theta.i[j] + h
      V.j <- SSN:::makeCovMat(
        theta.j, 
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
      ep[[j]] <- (V.j - V)/h
    }
    
    for(j in 1:np){
      for(k in j:np){
        I_REML <- P %*% ep[[j]] %*% P %*% ep[[k]]
        em[k, j] <- em[j, k] <- 1/2 * sum(diag(I_REML))
      }
    }
    totFIM <- em + old.information
    FIM[i] <- log(det(totFIM))
    
  }
  
  if(any(is.infinite(FIM))){
    return(-1e9)
  } else {
    return(mean(FIM, na.rm = TRUE))
  }
  
}
