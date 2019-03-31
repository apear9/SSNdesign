#' A convenience function that wraps \code{extractExpectedFisherInformation}
#' 
#' @description 
#' 
#' This function returns the expected standard errors on the covariance parameters for a \code{glmssn} object.
#' 
#' @param glmssn A fitted \code{glmssn} object.
#' @param log.scale Whether the standard errors should be given on the log scale. Defaults to FALSE.
#' @param ... Optional. Additional arguments to \code{extractExpectedFisherInformation}. 
#' @return A vector of expected standard errors for the \code{glmssn} object's covariance parameters.
#' @export
getExpectedSE <- function(glmssn, log.scale = FALSE, ...){
  fish <- extractExpectedFisherInformation(glmssn, ...)
  ster <- sqrt(diag(solve(fish)))
  if(log.scale){
    parm <- glmssn$estimates$theta
    parm <- parm[1:length(parm)]
    # Delta method transformation
    ster.log <- (1/parm) * ster 
    return(ster.log)
  } else {
    return(ster)
  }
}