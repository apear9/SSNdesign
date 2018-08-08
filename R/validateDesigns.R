#' Validate designs on SpatialStreamNetworks
#' 
#' @description 
#' 
#' This function performs validation on
#' 
#' @usage 
#' 
#' validateDesigns(designs.list, glmssn, type, dist = rmse, n.sims = 500)
#' 
#' @param full.ssn An object of class SpatialStreamNetwork from which the designs were originally selected.
#' @param designs.list A SpatialStreamNetwork objects containing the design points to be used for valdiation.
#' @param glmssn A fitted glmssn object 
#' @param type A string indicating what the design(s) should be tested on. This argument must be one of "fixed.effects", "covariance.parms" or "prediction".
#' @param dist A loss function comparing the simulated and predicted values at each step in the validation process. See Details for more information.
#' @param n.sims The number of times values should be simulated for the validation procedure. Defaults to 500.
#' @return A list of validation results. See Details for more information.
#' 
#' @details
#' 
#' The functions given to the argument dist should be functions of the form \code{function(x, y)}. They must at least accept two arguments. If the ordering of the arguments is important, then note that the simulated values are passed to the x argument and the predicted values to the y argument. 
#' 
#' The list of validation results returned by this function contains a matrix of dimensions m x n.sims, where m is equal to length(designs.list). The elements of this matrix represent a single validation result for one design at one simulation. It also contains a numeric vector of m elements, which are obtained by taking the rowmeans of the aforementioned matrix.
#' 
#' @export
validateDesigns <- function(full.ssn, designs.list, glmssn, type, dist.type = rmse, n.sims = 500){
  
  # Check that type argument is not longer than one
  if(length(type) > 1){
    warning("The argument type is longer than one, so only the first element of this vector will be used.")
    type <- type[1]
  }
  
  # Check that type argument is valid
  if(!(type %in% c("fixed.effects", "covariance.parms", "prediction"))){
    stop("The argument type must be one of fixed.effects, covariance.parms or prediction")
  }
  
  # Check what case to use
  results <- list()
  if(type == "fixed.effects"){
    results <- design.validation.fixed.effects.estimation(full.ssn = full.ssn, designs.list = designs.list, glmssn = glmssn, dist.type = dist.type, n.sims = n.sims)
  }
  if(type == "covariance.parms"){
    results <- design.validation.covariance.parms.estimation(full.ssn = full.ssn, designs.list = designs.list, glmssn = glmssn, dist.type = dist.type, n.sims = n.sims)
  }
  if(type == "prediction"){
    results <- design.validation.predictions(full.ssn = full.ssn, designs.list = designs.list, glmssn = glmssn, dist.type = dist.type, n.sims = n.sims)
  }
  
  # RETURN
  return(results)
  
}

