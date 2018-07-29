rmse <- function(x, y){
  
  diff <- x - y
  sqrt(mean(diff^2))
  
}