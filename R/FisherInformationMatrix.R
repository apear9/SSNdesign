#' A utility function for covariance parameter estimation
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
FisherInformationMatrix <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Rely on ordering of output of glmssn to retrieve correct rows in the design matrix
  ind <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  
  # Get design matrix -- NEEDS COLUMN OF ONES???
  X <- glmssn$sampinfo$X[ind, ]
  Xt <- t(X)
  cds <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind, ]
  colnames(cds) <- c("x", "y")
  
  ## Get distance matrices, etc.
  
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  
  ## Simulate covariance parameters
  
  cvp.cols <- length(glmssn$estimates$theta)
  cvp <- matrix(nrow = n.draws, ncol = cvp.cols)
  
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
  
  FIM <- vector("numeric", n.draws)
  h <- 1e-3
  for(i in 1:n.draws){
    
    # Get covariance matrix on the data
    
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
      ep[[j]] <- (V.j - V)/h
    }
    
    for(j in 1:np){
      for(k in j:np){
         I_REML <- P %*% ep[[j]] %*% P %*% ep[[k]]
         em[j, k] <- 1/2 * sum(diag(I_REML))
         em[k, j] <- em[j, k]
      }
    }
    print(em)
    iFIM.i <- solve(em)
    FIM[i] <- det(iFIM.i)
    
  }
  
  return(mean(FIM, na.rm = TRUE))
  
}
