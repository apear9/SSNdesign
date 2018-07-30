#' Extract covariates on stream edges to observed and/or predicted sites
#' 
#'@description 
#'
#' This function is designed to get covariates associated with the data slot of the SSN obejct and join them to the point.data slot of the observed and/or predicted sites slot.
#' 
#'@usage 
#' 
#'\code{extractStreamEdgeCovariates(ssn)}
#' 
#'@param ssn an object of class SpatialStreamNetwork
#'@param columns a vector of column names that should be extracted from the data slot of the SpatialStreamNetwork object
#'@return An object of class SpatialStreamNetwork. 
#'
#'@details
#'
#'This function can be useful when the stream edges contain data that you want to associate with the sites for modelling or prediction purposes. An example is watershed attributes for the stream segments.
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
