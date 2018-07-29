#' Overwrite values of the covariance parameters in a fitted glmssn object.
#' 
#' @description 
#' 
#' It is not recommended to manually assign values for the covariance parameters in a glmssn object, since the best estimates are already provided. However, it may be necessary for data simulation purposes.
#' 
#' @param glmssn An object of class glmssn.
#' @param new.theta A numeric vector containing the new values of the covariance parameters.
#' @return An obejct of class glmssn with the old values of theta replaced by the new ones.
#' 
#' @details 
#' 
#' Note that the length of new.theta must be equal to the length of \code{glmssn$estimates$theta} and the covariance parameter values must be provided in the same order.
#' 
#' @export
overwriteThetas <- function(glmssn, new.theta){
  n.old <- length(glmssn$estimates$theta)
  n.new <- length(new.theta)
  if(n.old != n.new){
    err <- paste("The new vector of covariance parameters must be of the same length as the old vector.\nCurrent state\nOld:", n.old, "\nNew:", n.new)
    stop(err)
  }
  glmssn$estimates$theta <- new.theta
  return(glmssn)
}