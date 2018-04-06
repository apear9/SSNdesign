#' Place a point at some distance along a stream edge
#'
#'@description 
#'   
#'\code{locatePointOnEdge()} returns the coordinates of a point along a line, as well as the distance from the first vertex to the point along the line.
#'
#'@usage 
#'
#'This function is not intended for direct use. It is called inside \code{\link{generateSites}} to ensure that the placement of simulated sites in geographical space is correct. 
#'   
#'@param ssn an object of class SpatialStreamNetwork
#'@param rid the unique identifier for an edge in a SpatialStreamNetwork
#'@param ratio the proportion of the total length of the stream edge where the point should be located
#'@return a vector containing the xy coordinates of the point, and its distance from the first point on the stream edge
#'
#'@export
locatePointOnEdge <- function(ssn, rid, ratio){
  # get line segments according to criteria
  ssncoords <- ssn@lines[[rid + 1]]@Lines[[1]]@coords
  nsegments <- nrow(ssncoords) - 1
  # get lengths of every segment
  seglengths <- rep(NA, nsegments)
  for(i in 1:nsegments){
    seglengths[i] <- sqrt(sum((ssncoords[i, ] - ssncoords[i + 1, ])^2))
  }
  sumseglengths <- cumsum(seglengths)
  edgelength <- max(sumseglengths)
  sumsegratios <- sumseglengths/edgelength
  # get coord of last vertex before distfromstart
  ids <- 2:(nsegments+1)
  last <- suppressWarnings(max(ids[sumsegratios <= ratio]))
  if(is.infinite(last)){
    last <- 1
  }
  # get coord of last vertex before and after distfromend
  lastvertex <- ssncoords[last, ]
  nextvertex <- ssncoords[last + 1, ]
  # calc as proportion of total edge length
#  lastprop1 <- sumseglengths[last]/edgelength
  loclength <- ratio * edgelength
  if(last - 1 != 0){
    scllength <- (loclength - sumseglengths[last - 1])/seglengths[last]
  } else {
    scllength <- loclength/seglengths[last]
  }
  
  # recalculate new proportion, being the proportion of the enclosed segment along which we will place the sampling point
#  if(!is.infinite(last))
#  newprop <- abs(ratio - lastprop1)/(seglengths[last])
    # calculate coordinate on the constrained edge
  newvertex <- lastvertex + scllength * (nextvertex - lastvertex)
  # return everything.
  return(c(newvertex, loclength))
}