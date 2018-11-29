#' Extract covariates on stream edges to observed and/or predicted sites
#' 
#'@description 
#'
#' This function transfers specified variables recorded on the edges of a SpatialStreamNetwork object to the observed (and possibly prediction) sites in that SpatialStreamNetwork object. 
#' 
#'@usage 
#' 
#'\code{extractStreamEdgeCovariates(ssn, columns)}
#' 
#'@param ssn an object of class SpatialStreamNetwork
#'@param columns a vector of column names that should be extracted from the data slot of the SpatialStreamNetwork object
#'@return An object of class SpatialStreamNetwork. 
#'
#'@examples
#'
#'set.seed(1)
#'
#'# Create an SSN
#'s <- createSSN(10, binomialDesign(10), path = paste(tempdir(), "s.ssn", sep = "/"), importToR = TRUE)
#'
#'# Extract stream edge covariates
#'s <- extractStreamEdgeCovariates(ssn, "Length")
#'
#'@export
extractStreamEdgeCovariates <- function(ssn, columns){
  
  # Number of columns
  nc <- length(columns)
  
  # Check inputs
  if(nc < 1){
    stop("You must specify at least one column to extract data from.")
  }
  
  # Check whether prediction sites exist
  preds.exist <- length(ssn@predpoints@SSNPoints) > 0
  # Add columns to observed +/- pred points and join
  obs.unique.rids <- unique(ssn@obspoints@SSNPoints[[1]]@point.data$rid)
  obs.unique.rids.ind <- match(obs.unique.rids, ssn@data$rid)
  # Get a list of values in each column
  col.list <- vector("list", nc)
  for(i in 1:nc){
    col.list[[i]] <- ssn@data[obs.unique.rids.ind, columns[i]]
  }
  names(col.list) <- columns
  n.unique.rids <- length(obs.unique.rids)
  for(i in 1:n.unique.rids){
    for(j in 1:nc){
      ind <- which(ssn@obspoints@SSNPoints[[1]]@point.data$rid == obs.unique.rids[i])
      ssn@obspoints@SSNPoints[[1]]@point.data[ind, columns[[j]]] <- col.list[[j]][i]
    }
  }
  if(preds.exist){
    pred.unique.rids <- unique(ssn@predpoints@SSNPoints[[1]]@point.data$rid)
    pred.unique.rids.ind <- match(pred.unique.rids, ssn@data$rid)
    col.list <- vector("list", nc)
    for(i in 1:nc){
      col.list[[i]] <- ssn@data[pred.unique.rids.ind, columns[i]]
    }
    names(col.list) <- columns
    n.unique.rids <- length(pred.unique.rids)
    for(i in 1:n.unique.rids){
      for(j in 1:nc){
        ind <- which(ssn@predpoints@SSNPoints[[1]]@point.data$rid == pred.unique.rids[i])
        ssn@predpoints@SSNPoints[[1]]@point.data[ind, columns[[j]]] <- col.list[[j]][i]
      }
    }
  }
  
  # Return
  return(ssn)
  
}
