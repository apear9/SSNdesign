#' @inherit DOptimality
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
