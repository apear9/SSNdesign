#' Wrapper function for SimulateOnSSN
#' 
#' @description This function was developed as a workaround for an error that occurs in \code{SSN::SimulateOnSSN} when a \code{SpatialStreamNetwork} object has temporal replicates on its sites.
#' 
#' @param ... Arguments for \code{SSN::SimulateOnSSN}.
#' @return A \code{SpatiaLStreamNetwork} object like the one returned in the \code{ssn.object} element of the list output from \code{SSN::SimulateOnSSN}. 
#' 
#' @export
SimulateOnSSN2 <- function(...){
  
  # Simulate data as usual (but keep only the ssn object)
  ssn <- SimulateOnSSN(...)$ssn.object
  
  # Deal with the obs locID problem
  obs <- getSSNdata.frame(ssn)
  obs$locID <- as.character(obs$locID) # force it to be character
  obs$locID <- gsub("o", "", obs$locID) # remove the o
  obs$locID <- factor(
    as.numeric(obs$locID)
  ) # force it to be numeric, then factor
  ssn@obspoints@SSNPoints[[1]]@point.data <- obs
  
  # Same for preds if any are present
  if(length(ssn@predpoints@SSNPoints) > 0){
    prd <- getSSNdata.frame(ssn, ssn@predpoints@ID[[1]])
    prd$locID <- as.character(prd$locID)
    prd$locID <- gsub("p", "", prd$locID)
    prd$locID <- factor(
      as.numeric(prd$locID)
    ) 
    ssn@predpoints@SSNPoints[[1]]@point.data <- prd
  }
  
  ssn
  
}
