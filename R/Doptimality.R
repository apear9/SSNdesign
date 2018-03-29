#' A utility function for fixed effects parameter estimation with known covariance parameters
#' 
#' \code{DOptimality()}
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param glmssn An object of class glm.ssn
#' @param design.points A numeric vector of pids for sites in a proposed design
#' @param prior.parameters A list of functions parameterised in terms of n.draws that simulate from probability distributions
#' @param extra.arguments A list of extra arguments that the utility may rely on
#' @return A numeric representing the D-optimal utility of a design estimated by Monte-Carlo integration
#' 
#' @export
DOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # OBSOLETE CODE
  # ## Get distance matrices, etc.
  # 
  # distance.matrix <- important.matrices$Distance
  # connectivity.matrix <- important.matrices$Connectivity
  # weights.vector <- rep(1, nrow(connectivity.matrix))
  # conn.weight.matrix <- connectivity.matrix 
  # 
  # ## Get other arguments
  # 
  # cvp <- vector("list", length(prior.parameters)) 
  # n <- nrow(distance.matrix)
  # X <- model.matrix(object = formula, data = design)
  # Xt <- t(X)
  
  ## Get data for design points
  
  # print(design.points)
  
  # Rely on ordering of output of glmssn to retrieve correct rows in the design matrix
  ind <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  
  # print(ind)
  
  X <- glmssn$sampinfo$X[ind, ]
  
  # print(X)
  
  Xt <- t(X)
  cds <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind, ]
  # print(class(cds))
  # 
  colnames(cds) <- c("x", "y")
  
  ## Get afv for design points
  # Not needed in this version of the function
  # afvColumn <- glmssn$args$addfunccol
  # design.afv <- ssn@obspoints@SSNPoints[[1]]@point.data[ind, afvColumn]
  # 
  ## Get distance matrices, etc.
  
  #dist.junc.obs <- SSN:::getStreamDistMatInt(ssn, design.points, "obs")
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  
  # print(mat)
  
  # This code is OBSOLETE
  # distance.matrix <- important.matrices$Distance
  # connectivity.matrix <- important.matrices$Connectivity
  # weights.vector <- rep(1, nrow(connectivity.matrix))
  # conn.weight.matrix <- connectivity.matrix 
  # 
  ## Get simulated covariance parameters
  
  cvp.cols <- length(glmssn$estimates$theta)
  cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  
  # OBSOLETE CODE
  # X <- model.matrix(object = formula, data = design)
  # Xt <- t(X)
  
  ## Get other model parameters
  
  td <- glmssn$args$useTailDownWeight
  cm <- glmssn$args$CorModels
  un <- glmssn$args$use.nugget
  ua <- glmssn$args$use.anisotropy
  re <- glmssn$sampInfo$REs
  
  ## Simulate covparms from priors
  
  for(i in 1:cvp.cols){
    
    cvp[, i] <- prior.parameters[[i]](n.draws)
    
  }
  
  ## Perform MC simulations
  
  D <- vector("numeric", n.draws)
  
  for(i in 1:n.draws){
    
    theta.i <- cvp[i, ]
    # print(theta.i)
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
    
    ##	if(class(qrV) != "qr") {browser()}
    covbi <- t(X) %*% solve(V) %*% X
    covb <- solve(covbi)
    D[i] <- -log(det(covb))
    
  }
  
  return(mean(D, na.rm = TRUE))
  
}
