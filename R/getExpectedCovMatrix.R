#' A convenience function that wraps \code{extractExpectedFisherInformation}
#' 
#' @description 
#' 
#' This function returns the expected covariance matrix for the covariance parameters in a \code{glmssn} object. This function is similar to \code{\link{getExpectedSE}}.
#' 
#' @param glmssn A fitted \code{glmssn} object.
#' @param log.scale Whether the covariance matrix should be on the log-scale. Defaults to FALSE.
#' @param ... Optional. Additional arguments to \code{extractExpectedFisherInformation}.
#' @return A vector of expected standard errors for the glmssn object's covariance parameters.
#' @export
getExpectedCovMatrix <- function(glmssn, log.scale = FALSE, ...){
  fish <- extractExpectedFisherInformation(glmssn, ...)
  cov <- solve(fish)
  if(log.scale){
    parm <- glmssn$estimates$theta
    parm <- parm[1:length(parm)]
    I <- nrow(cov)
    J <- ncol(cov)
    cov.log <- matrix(0, I, J)
    for(i in 1:I){
      for(j in 1:J){
        # Delta method transformation to log-scale
        cov.log[i, j] <- 1/(parm[i] * parm[j]) * cov[i, j]
      }
    }
    return(cov.log)
  } else {
    return(cov)
  }
}