#' Transform to standard normal
#' 
#' @description This function takes the vector x and transforms it to be approximately standard normal. It is meant to be used with the function \code{\link{transformSSNVars}}.
#' 
#'@param x A numeric vector of observations
#'@param na.rm Whether missing observations should be ignored. Defaults to \code{TRUE}. 
#'@return A numeric vector of size \code{length(x)}.
#'
#'@examples
#'
#'x <- rnorm(100, 4, 2)
#'x.std <- stdsGs(x)
#'par(mfrow = c(1, 2))
#'plot(density(x), main = "Unstandardised")
#'plot(density(x.std), main = "Standardised")
#'
#'@export
stdsGs <- function(x, na.rm = T){
  (x - mean(x, na.rm = na.rm))/sd(x, na.rm = na.rm)
}
