#' A function to solve adaptive design problems.
#' 
#' @description 
#' 
#' This function is a wrapper for \code{\link{optimiseSSNDesign}}, \code{\link{splitSSNSites}} and \code{\link{spliceSSNSites}}. It uses these functions in a way that automates the adaptive design process.
#' 
#' The adaptive design process used in this function is called myopic design (). 
#' 
#' @usage 
#' 
#' doAdaptiveDesign(ssn, new.path, n.points.per.step, legacy.sites.step, fixed.once.chosen = TRUE, utility.function, prior.parameters, ...)
#' 
#' @param ssn An object of class SpatialStreamNetwork. This network must have proposed sites for several time steps.
#' @param new.path A path to a folder where new SpatialStreamNetwork objects and intermediate results will be stored. Note that this folder should be a .ssn directory. See Details for more information about how files will be stored in this folder. 
#' @param glmssn An object of class \code{glmssn}. This must be the 'true' model that will be used in all adaptive design steps.
#' @param n.points.per.step A list of vectors, which all resemble the \code{n.points} argument required by \code{\link{optimiseSSNDesign}}.
#' @param legacy.sites.step A list of vectors, which all resemble the \code{legacy.sites} argument in \code{\link{optimiseSSNDesign}}.
#' @param fixed.once.chosen A logical. Defaults to \code{TRUE}.
#' @param utility.function The utility function which should be used to select the design at each step. Only one utility function can be specified. 
#' @param prior.parameters A list of 
#' @param ... Additional arguments to \code{\link{optimiseSSNDesign}}.
#' @return A list of SSNdesign lists, which are the outputs of \code{\link{optimiseSSNDesign}}. There will be ... Any files will also remain in the folder specified by \code{new.path}.  
#'
#'  @details 
#'  
#'  Something which is no doubt very important.
#'  
#'  @examples 
#'  
#'  \dontrun{#CODE}
#'  
#'  @export
doAdaptiveDesign <- function(ssn, new.path, n.points.per.step, legacy.sites.step, fixed.once.chosen = TRUE, utility.function, prior.parameters, ...){
  
  NULL
  
}