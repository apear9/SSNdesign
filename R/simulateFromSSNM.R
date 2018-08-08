#' Simulate data, with error, from a glmssn object
#' 
#' @description 
#' 
#' This function simulates data on a SpatialStreamNetwork based on information from a fitted model object. It essentially acts as a wrapper for SimulateOnSSN. 
#' 
#' @usage 
#' 
#' \code{simulateFromSSNM(ssn, glmssn, fixed.effects = NULL, covariance.type = NULL, covariance.parms = NULL)}
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param glmssn An object of class glmssn
#' @param fixed.effects A numeric vector specifying the values of the fixed effects. Leave as NULL unless you specifically intend to override the information about the fixed effects you get from the glmssn.
#' @param covariance.type A character vector specifying what kind of covariance structure should be used to simulate the data. Leave as NULL unless you specifically intend to overwrite the covariance components from the glmssn. 
#' @param covariance.parms A numeric vector specifiying the values of the covariance parameters. Leave as NULL unless you specifically intend to override the informaiton about the covariance parameters you get from the glmssn. 
#' @return An object of class SSNM, which is a modified version of the one that was given as the argument ssn.
#' 
#' @details 
#' 
#' It is only possible to use this function to simulate data on a SpatialStreamNetwork object when a fitted model exists.
#' 
#' @export
simulateFromSSNM <- function(ssn, glmssn, fixed.effects = NULL, covariance.type = NULL, covariance.parms = NULL){
  
  # Check whether we need to overwrite the fixed effects from the model
  if(!is.null(fixed.effects)){
    if(!is.numeric(fixed.effects)){
      stop("The argument fixed.effects must be either NULL or numeric.")
    }
    glmssn$estimates$betahat <- fixed.effects
  }
  
  # Check whether we need to overwrite the covariance function as a whole
  if(!is.null(covariance.type)){
    if(!is.character(covariance.type)){
      stop("The argument covariance.type must be either NULL or character.")
    }
    glmssn$args$CorModels <- covariance.type
  }
  
  # Check whether we need to do the same for the covariance parameters
  if(!is.null(covariance.parms)){
    if(!is.numeric(covariance.parms)){
      stop("The argument fixed.effects must be either NULL or numeric.")
    }
    glmssn$estimates$theta <- covariance.parms
  }
  
  # Rewrite model formula as one sided formula
  form <- as.character(glmssn$args$formula)
  form <- paste(form[1], form[3])
  form <- as.formula(form)
  
  # Get dataframes from the SSN
  o.df <- getSSNdata.frame(ssn)
  if(anyPreds(ssn)){
    p.df <- getSSNdata.frame(ssn, "preds")
    preds <- "preds"
  } else {
    p.df <- NULL
    preds <- NULL
  }
  
  # Simulate data to SpatialStreamNetwork
  ga <- glmssn$args #arguments
  ge <- glmssn$estimates # estimates
  
  ssn.new <- SimulateOnSSN(
    ssn,
    o.df,
    p.df,
    preds,
    form,
    ge$betahat,
    ga$CorModels,
    ga$use.nugget,
    ga$use.anisotropy,
    ge$theta,
    ga$addfunccol,
    ga$useTailDownWeight,
    ga$family,
    FALSE
  )$ssn.object
  
  ssn.new@obspoints@SSNPoints[[1]]@point.data$locID <- as.character(
    ssn.new@obspoints@SSNPoints[[1]]@point.data$locID
  )
  ssn.new@obspoints@SSNPoints[[1]]@point.data$locID <- gsub(
    "o",
    "",
    ssn.new@obspoints@SSNPoints[[1]]@point.data$locID
  )
  ssn.new@obspoints@SSNPoints[[1]]@point.data$locID <- anum(
    ssn.new@obspoints@SSNPoints[[1]]@point.data$locID
  )
  if(!is.null(preds)){
    ssn.new@predpoints@SSNPoints[[1]]@point.data$locID <- as.character(
      ssn.new@predpoints@SSNPoints[[1]]@point.data$locID
    )
    ssn.new@predpoints@SSNPoints[[1]]@point.data$locID <- gsub(
      "p",
      "",
      ssn.new@predpoints@SSNPoints[[1]]@point.data$locID
    )
    ssn.new@predpoints@SSNPoints[[1]]@point.data$locID <- anum(
      ssn.new@predpoints@SSNPoints[[1]]@point.data$locID
    )
  }
  
  # Return result
  return(ssn.new)
  
}