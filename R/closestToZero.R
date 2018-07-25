closestToZero <- function(d,x){
  # d is an element in a list of distances
  # x is a row of a hydrological distance matrix
  # X MUST BE NAMED
  if(is.null(names(x))){
    stop("x must be a named vector, representing a row from a hydrological distance matrix")
  }
  x <- abs(x - d)
  cTZ <- which(x == min(x))
  return(names(cTZ))
}
