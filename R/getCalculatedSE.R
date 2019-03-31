#' A convenience function that extracts the standard errors for the covariance parameters from a \code{glmssn} object
#' 
#' @description 
#' 
#' This function returns the expected standard errors on the covariance parameters for a given set of parameter values.
#' 
#' @param glmssn A fitted \code{glmssn} object.
#' @param log.scale Whether the standard errors should be on the log-scale. Defaults to \code{FALSE}.
#' @return A vector of expected standard errors for the glmssn object's covariance parameters.
#' @export
getCalculatedSE <- function(glmssn, log.scale = FALSE){
  
  fish.exists <- !is.null(glmssn$optimOutput$hessian)
  if(fish.exists){
    fish <- glmssn$optimOutput$hessian
  } else {
    stop("Fisher Information Matrix does not exist. Suggestion: use getExpectedSE() instead.")
  }
  
  ster.log <- sqrt(diag(solve(fish)))
  if(log.scale){
    return(ster.log)
  } else {
    parm <- glmssn$estimates$theta
    parm <- parm[1:length(parm)]
    ster <- parm * ster.log # undo delta method transformation
    return(ster)
  }
  
}