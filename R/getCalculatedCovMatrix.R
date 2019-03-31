#' A convenience function that extracts the standard errors for the covariance parameters from a \code{glmssn} object
#' 
#' @description 
#' 
#' This function returns the estimated covariance matrix on the covariance parameters from a \code{glmssn} object. Importantly it can be set the return the estimated covariance matrix on its native scale, as opposed to the log scale. 
#' 
#' @param glmssn A fitted \code{glmssn} object.
#' @param log.scale Whether the matrix should be on the log-scale. Defaults to \code{FALSE}.
#' @return A covariance matrix for the \code{glmssn} object's covariance parameters.
#' @export
getCalculatedCovMatrix <- function(glmssn, log.scale = FALSE){
  
  fish.exists <- !is.null(glmssn$optimOutput$hessian)
  if(fish.exists){
    fish <- glmssn$optimOutput$hessian
  } else {
    stop("Fisher Information Matrix does not exist. Suggestion: use getExpectedCovMatrix() instead.")
  }
  
  cov.log <- solve(fish)
  if(log.scale){
    return(cov.log)
  } else {
    parm <- glmssn$estimates$theta
    parm <- parm[1:length(parm)]
    I <- nrow(cov.log)
    J <- ncol(cov.log)
    cov <- matrix(0, I, J)
    for(i in 1:I){
      for(j in 1:J){
        # Undo delta method transformation to log-scale
        cov[i, j] <- parm[i] * parm[j] * cov.log[i, j]
      }
    }
    return(cov)
  }
  
}