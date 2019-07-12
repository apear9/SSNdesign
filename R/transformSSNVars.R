#' A function to conveniently transform multiple variables inside the observed and predicted points' data slots in a SpatialStreamNetwork object
#' 
#' @description 
#' 
#' The process of modifying the data.frame objects in the observed and predicted sites of a SpatialStreamNetwork can be complicated. This function simplifies this process.
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param new.ssn.path The directory in which a new .ssn directory should be written out. Ignored if write.out = FALSE.
#' @param obs.vars A character vector containing the names of the columns in the observed sites point.data slot which should be transformed. If this is a named vector, then new columns with those names will be added. Otherwise, the old columns will be overwritten.
#' @param preds.vars A character vector containing the names of the columns in the prediction sites point.data slot which should be transformed. If this is a named vector, then new columns with those names will be added. Otherwise, the old columns will be overwritten.
#' @param func The function which should be used to transform the columns.
#' @param write.out A logical indicating whether the .ssn folder in \code{ssn@path} should be updated.
#' @param ... Additional arguments to func. 
#' @return An object of class SpatialStreamNetwork
#' 
#' @details 
#' 
#' Note, if write.out = TRUE then the SSN will be re-imported before returning. Note also that write.out = TRUE will cause the data inside the .ssn directory to be overwritten. 
#' 
#' @export
transformSSNVars <- function(ssn, new.ssn.path, obs.vars, preds.vars, func, write.out = FALSE, ...){
  
  # Basic argument checking
  if(!is.function(func)){
    stop("The argument func must be a function.")
  }
  if(!is.logical(write.out)){
    stop("The argument write.out must be a logical. Defaults to FALSE.")
  }
  if(missing(obs.vars) & missing(preds.vars)){
    stop("At least one of obs.vars or preds.vars must be specified.")
  }

  # Compute values for the observed sites
  if(!missing(obs.vars)){
    obs <- getSSNdata.frame(ssn)
    if(is.null(names(obs.vars))){
      obs[,obs.vars] <- apply(as.matrix(obs[,obs.vars]), 2, func, ...)
    } else {
      to.add <- apply(as.matrix(obs[,obs.vars]), 2, func, ...)
      colnames(to.add) <- names(obs.vars)[1:length(obs.vars)]
      obs <- cbind(obs, to.add)
    }
    ssn <- putSSNdata.frame(obs, ssn)
  }
  
  # Same for preds
  if(!missing(preds.vars)){
    prd <- getSSNdata.frame(ssn, "preds")
    if(is.null(names(preds.vars))){
      prd[,preds.vars] <- apply(as.matrix(prd[,preds.vars]), 2, func, ...)
    } else {
      to.add <- apply(as.matrix(prd[,preds.vars]), 2, func, ...)
      colnames(to.add) <- names(preds.vars)[1:length(preds.vars)]
      prd <- cbind(prd, to.add)
    }
    ssn <- putSSNdata.frame(prd, ssn, "preds")
  }
  
  # Write out to file if asked
  if(write.out){
    writeSSN(ssn, new.ssn.path)
    if(anyPreds(ssn)){
      preds <- "preds"
      ssn <- importSSN(new.ssn.path, preds)
    } else {
      ssn <- importSSN(new.ssn.path)
    }
  }
  
  # Return ssn object
  return(ssn)
  
}

