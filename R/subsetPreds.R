#' Subset the prediction sites in a SpatialStreamNetwork
#' 
#' @description 
#' 
#' This function works like subsetSSN, but instead of focussing on the observed sites, this allows the prediction sites alone to be subset according to some logical criteria. It is particularly useful when attempting to keep only prediction sites with certain PIDs.
#' 
#' @usage
#' 
#' \code{subsetPreds(ssn, subset)}
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param subset A logical condition used to subset the prediction sites.
#' @return A SpatialStreamNetwork. The new preds shapefile overwrites the old one in \code{ssn@path}.
#' 
#' @export
subsetPreds <- function(ssn, subset){
  
  pred.len <- length(ssn@predpoints@SSNPoints)
  
  ## Meat from subsetSSN from SSN
  if (pred.len > 0) {
    for (i in 1:pred.len) {
      pred.name <- ssn@predpoints@ID[[i]]
      ind.preds <- eval(substitute(subset), ssn@predpoints@SSNPoints[[i]]@point.data)
      ind.na <- is.na(ind.preds)
      ind.preds[ind.na] <- FALSE
      rm(ind.na)
      
      if(sum(ind.preds)==1)
        coords<- as.matrix(t(ssn@predpoints@SSNPoints[[i]]@point.coords[ind.preds, ]))
      else coords<- ssn@predpoints@SSNPoints[[i]]@point.coords[ind.preds,]
      proj4string <- ssn@proj4string
      data.tmp <- ssn@predpoints@SSNPoints[[i]]@point.data[ind.preds,]
      
      ind.xy <- names(data.tmp) == "coords_x1" | names(data.tmp) == "coords_x2"
      if (sum(ind.xy) > 0) {
        data.tmp <- data.tmp[,!ind.xy]}
      
      preds.sub <- SpatialPointsDataFrame(coords = coords, data = data.tmp, proj4string = proj4string)
      writeOGR(preds.sub, paste0(ssn@path,"/",pred.name, ".shp"), pred.name, "ESRI Shapefile", overwrite_layer = TRUE)
      SSN:::write.dbf.SSN(preds.sub@data, paste0(ssn@path,"/",pred.name), max_nchar = 30)
      rm(coords, proj4string, data.tmp, preds.sub, ind.preds, ind.xy)
    }
  } else {
    warning("No prediction points found! No changes have been made to the SSN object.")
    return(ssn)
  }
  
  ssn <- importSSN(ssn@path, predpts = pred.name[1], o.write = F)
  return(ssn)
  
}
