#' Overwrite values of the fixed effects parameters in a fitted glmssn object.
#' 
#' @description 
#' 
#' It is not recommended to manually assign values for the fixed effects parameters in a glmssn object, since the best estimates are already provided. However, it may be necessary for data simulation purposes.
#' 
#' @param glmssn An object of class glmssn.
#' @param new.beta A numeric vector containing the new values of the fixed effects parameters.
#' @return An obejct of class glmssn with the old values of theta replaced by the new ones.
#' 
#' @details 
#' 
#' Note that the length of new.beta must be equal to the length of \code{glmssn$estimates$betahat} and the fixed effects parameter values must be provided in the same order.
#' 
#' @export
overwriteBeta <- function(glmssn, new.beta){
  n.old <- length(glmssn$estimates$betahat)
  n.new <- length(new.beta)
  if(n.old != n.new){
    err <- paste("The new vector of fixed effects parameters must be of the same length as the old vector.\nCurrent state\nOld:", n.old, "\nNew:", n.new)
    stop(err)
  }
  glmssn$estimates$betahat <- new.beta
  return(glmssn)
}