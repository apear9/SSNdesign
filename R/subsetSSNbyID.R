#' Subset the observed and/or prediction sites by ID
#' 
#' @description 
#' 
#' This is a function to subset the observed and/or prediction sites on a SpatialStreamNetwork based on pids or locIDs.
#' 
#' @usage 
#' 
#' subsetByID(ssn, new.ssn.path, obs, preds, locID = FALSE)
#' 
#' @param ssn An object of class SpatialStreamNetwork.
#' @param new.ssn.path A path to a new .ssn directory where the outputs of this function should be stored. 
#' @param obs A numeric or character vector containing the pids or locIDs of the observed sites that should be kept in the subset.
#' @param preds A numeric or character vector containing the pids or locIDs of the prediction sites that should be kept in the subset. This argument can be skipped.
#' @param locID A logical indicating whether the obs and preds vectors contain pids or locIDs. Defaults to FALSE.
#' @return An object of class SpatialStreamNetwork containing only the specified subset of sites.
#' 
#' @details 
#' 
#' Note, the SpatialStreamNetwork that is returned will have been re-imported from the new.ssn.path directory after the subset operations are complete.
#' 
#' @export
subsetByID <- function(ssn, new.ssn.path, obs, preds, locID = FALSE){
  
  # Check inputs and clean
  if(missing(preds)){
    preds <- NULL
  }
  if(!is.logical(locID)){
    stop("The argument locID must be logical. Set TRUE to use locIDs for the subsets.")
  }
  
  # Split by case
  if(locID){
    result <- subsetByLOCID(ssn, new.ssn.path, obs, preds)
  } else {
    result <- subsetByPID(ssn, new.ssn.path, obs, preds)
  }
  
  # Re-import just to be sure that everything is coded correctly
  result <- suppressMessages(importSSN(new.ssn.path, preds))
  
  # Spit out result
  return(result)
  
}