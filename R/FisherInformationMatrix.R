FisherInformationMatrix <- function(formula, design, important.matrices, prior.parameters, n.draws, extra.arguments){
  
  distance.matrix <- important.matrices$Distance
  connectivity.matrix <- important.matrices$Connectivity
  weights.matrix <- matrix(data = 1, nrow = nrow(distance.matrix), ncol = ncol(distance.matrix))
  n <- nrow(distance.matrix)
  parms.sim <- lapply(prior.parameters, FUN = function(x){x(n.draws)})
  parms.sim.ev <- vector("list", n.draws)
  for(i in 1:n.draws){
    parms.sim.ev[[i]] <- c(
      parms.sim[[1]][i],
      parms.sim[[2]][i],
      parms.sim[[3]][i],
      parms.sim[[4]][i],
      parms.sim[[5]][i]
    )
  }
  parms.sim <- parms.sim.ev
  
  results <- lapply(
    parms.sim,
    FUN = function(x, formula, data_, distance.matrix, connectivity.matrix, weights.matrix, n, covariance.function){
      
      nk <- length(x)
      X <- model.matrix(formula, data_)
      Xt <-t(X)
      W <- covariance.function(x, distance.matrix, connectivity.matrix, weights.matrix, n)
      Wi <- solve(W)
      nk <- length(parms.sim)
      hv <- rep(0, nk)
      dWdk <- vector("list", nk)
      for(i in 1:nk){
        hv[i] <- 1e-6
        dWdk[[i]] <- (covariance.function(x + hv, distance.matrix, connectivity.matrix, weights.matrix, n) - W) / 1e-6
      }
      P <- Wi - Wi %*% X %*% solve(Xt %*% Wi %*% X) %*% Xt %*% Wi
      ij <- expand.grid(1:nk, 1:nk)
      FIM.ij <- apply(ij, 1, FUN = function(x){1/2 * trace( P %*% dWdk[[x[1]]] %*% P %*% dWdk[[x[2]]])})
      FIM <- matrix(nrow = nk, ncol = nk)
      for(i in 1:nk){
        
        for(j in 1:nk){
          
          FIM[i, j] <- FIM.ij[ij[, 1] == i, ij[, 2] == j]
          
        }
        
      }
      
      return(det(FIM))
      
    },
    formula = formula,
    data_ = design,
    distance.matrix = distance.matrix,
    connectivity.matrix = connectivity.matrix,
    weights.matrix = weights.matrix,
    n = n,
    covariance.function = covariance.function
  )
  
  results <- unlist(results)
  U_ <- mean(results)
  return(U_)
  
}