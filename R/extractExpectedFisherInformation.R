#' Extract the expected Fisher information matrix for the covariance parameters in a glmssn object
#' 
#' @description
#' 
#' This function takes a \code{glmssn} object and returns the expected Fisher information matrix on its covariance parameters.
#' 
#' @param glmssn An object of class \code{glmssn}.
#' @param covparms A vector of covariance parameter values. These are extracted from the glmssn object by default.
#' @param diff The step size used for finite differencing. 
#' @return The expected Fisher information matrix.
extractExpectedFisherInformation <- function(glmssn, covparms = glmssn$estimates$theta, diff = 1e-8){
  
  # Get information ready
  ssn <- glmssn$ssn.object
  cds <- ssn@obspoints@SSNPoints[[1]]@point.coords
  X <- glmssn$sampinfo$X
  Xt <- t(X)
  
  # Extract matrices needed for calculations
  d.m <- getStreamDistMatInOrder(ssn)
  d.j <- constructTotalMatrix(d.m)
  n.z <- d.j$net.zero
  d.j <- d.j$d.junc
  if(!is.null(glmssn$args$addfunccol)){
    addfunccol <- unlist(unname(getSSNdata.frame(ssn)[, glmssn$args$addfunccol]))
  }else{
    n.r <- nrow(getSSNdata.frame(ssn))
    addfunccol <- rep(1, n.r)
  }
  im <- getImportantMatrices.obs(d.j, addfunccol)
  the <- covparms
  the <- the[1:length(the)]
  cds <- ssn@obspoints@SSNPoints[[1]]@point.coords
  V <- SSN:::makeCovMat(
    the, 
    im$d, 
    im$a, 
    im$b, 
    im$w, 
    n.z, 
    cds[, 1], 
    cds[, 2], 
    cds[, 1], 
    cds[, 2], 
    glmssn$args$useTailDownWeight,
    glmssn$args$CorModels, 
    glmssn$args$use.nugget,
    glmssn$args$use.anisotropy,
    glmssn$sampinfo$REs
  )
  Vi <- solve(V)
  #	precompute the matrix P
  P <- Vi - Vi %*% X %*% solve(Xt %*% Vi %*% X) %*% Xt %*% Vi
  # estimate partial derivatives
  np <- length(the) # np = number of parameters
  ep <- vector("list", np) # empty vector of partial derivatives (hence e p)
  em <- matrix(0, nrow = np, ncol = np) # empty vector of matrices (hence e m)
  h <- diff
  for(j in 1:np){
    theta.j <- the
    theta.j[j] <- the[j] + h
    V.j <- SSN:::makeCovMat(
      theta.j, 
      im$d, 
      im$a, 
      im$b, 
      im$w, 
      n.z, 
      cds[, 1], 
      cds[, 2], 
      cds[, 1], 
      cds[, 2], 
      glmssn$args$useTailDownWeight,
      glmssn$args$CorModels, 
      glmssn$args$use.nugget,
      glmssn$args$use.anisotropy,
      glmssn$sampinfo$REs
    )
    ep[[j]] <- (V.j - V)/h
  }
  
  for(j in 1:np){
    for(k in j:np){
      I_REML <- P %*% ep[[j]] %*% P %*% ep[[k]]
      em[k, j] <- em[j, k] <- 1/2 * sum(diag(I_REML))
    }
  }
  
  # Return expected Fisher information
  em
  
}