stdsGs <- function(x, na.rm = T){
  (x - mean(x, na.rm = na.rm))/sd(x, na.rm = na.rm)
}