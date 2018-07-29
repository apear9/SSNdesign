anyPreds <- function(x){
  if(class(x)[1] != "SpatialStreamNetwork"){
    stop("The class of x must be a SpatialStreamNetwork")
  }
  length(
    x@predpoints@SSNPoints
  ) > 0
}
