#' Generate paths for .ssn folders to be written to tempdir()
#' 
#' @description This function takes the name of a .ssn folder and returns the full path to that .ssn folder in R's temporary directory.
#' 
#' @param ssn.name The name of the ssn folder (including .ssn extension) as a string.
#' @return The path for the \code{ssn.name} folder in R's temporary directory.
#' 
#' @examples
#' 
#' new.path <- tempPath("s.ssn")
#' s <- createSSN(5, binomialDesign(10), path = new.path, importToR = TRUE)
#' 
#' @export
tempPath <- function(ssn.name){
  
  if(!is.character(ssn.name)) stop("ssn.name must be a string.")
  
  return(paste(tempdir(), ssn.name, sep = "/"))
  
}