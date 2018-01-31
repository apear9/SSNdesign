DOptimality <- function(formula, design, important.matrices, prior.parameters, n.draws, extra.arguments){
  
  distance.matrix <- important.matrices$Distance
  connectivity.matrix <- important.matrices$Connectivity
  weights.matrix <- matrix(data = 1, nrow = nrow(distance.matrix), ncol = ncol(distance.matrix))
  cvp <- extra.arguments$covparms
  n <- nrow(distance.matrix)
  uCov <- cvp[1] * connectivity.matrix * weights.matrix * exp(-3 * distance.matrix / cvp[2]) 
  dCov <- cvp[3] * exp(-3 * distance.matrix / cvp[4])
  W  <- uCov + dCov + cvp[5] * diag(n)
  Wi <- solve(W)
  X <- model.matrix(object = formula, data = design)
  Xt <- t(X)
  D <- solve(Xt %*% Wi %*% X)
  D <- det(D)
  return(D)
  
}
