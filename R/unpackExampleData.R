#' Copy package data to a new folder
#' 
#' @description The SSNdesign package contains several .ssn and .Rdata files. These are the files referred to in the package vignette. This function takes these files out of SSNdesign's folders and copies them to a user-specified location on the user's computer.
#' 
#' @param unpack.loc The path to some folder on the user's computer that the files should be copied to.
#' @return A logical indicating whether file transfer was successful.
#' 
#' @details Please note that the function will return \code{FALSE}, indicating unsuccessful file transfer, if the directory specified in \code{unpack.loc} already contains files with the same names as the example data files. 
#' 
#' @examples
#'
#'\dontrun{
#'
#'# To the documents folder:
#'unpackExampleData("~/")
#'}
#' 
#' @export
unpackExampleData <- function(unpack.loc){
  if(!dir.exists(unpack.loc)) stop("The specified directory does not exist. Check the input for mistakes.")
  current.folder <- system.file("extdata", package = "SSNdesign")
  packedup.files <- dir(current.folder, full.names = TRUE)
  file.copy(packedup.files, unpack.loc)
}