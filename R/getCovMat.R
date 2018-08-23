#' Extract out the covariance matrix for a SpatialStreamNetwork for an arbitrary covariance structure
#' 
#' @description 
#' 
#' It is sometimes useful to extract out the covariance matrix associated with a covariance structure to test whether it is viable for the data. This function simplifies the process of doing so.
#' 
#' @param ssn An object of class SpatialStreamNetwork
#' @param glmssn A fitted glmssn object. Optional if other arguments are provided.
#' @param addfunccol The name of the column in the point.data slot of the obspoints slot of a SpatialStreamNetwork object containing the additive function values used to weight the tail up component of a covariance function.
#' @param theta A vector of covariance parameters. Optional if glmssn is provided.
#' @param useTailDownWeight A logical indicating whether the taildown component of a covariance function should be weighted by the additive function value
#' @param CorModels The covariance components to use.
#' @param use.nugget A logical indicating whether a nugget effect should be included in the covariance function.
#' @param use.anisotropy NOT YET IMPLEMENTED IN SSN. SET TO FALSE IF REQUIRED.
#' @param REs A logical indicating whether random effects are to be used. 
#' @return A list containing the covariance matrix and its determinant.
#' 
#' @export
getCovMat <- function(ssn, glmssn, addfunccol, ...){
  # Get coords
  cds <- ssn@obspoints@SSNPoints[[1]]@point.coords
  x <- cds[,1]
  y <- cds[,2]
  # Get important matrices
  d.junc.list <- getStreamDistMatInOrder(ssn)
  total.matrix <- constructTotalMatrix(d.junc.list)
  n.zero <- total.matrix$net.zero
  d.junc <- total.matrix$d.junc
  # Construct covariance matrix
  if(is.missing(glmssn)){
    i.m <- getImportantMatrices.obs(d.junc, ssn@obspoints@SSNPoints[[1]]@point.data[, addfunccol])
    cov.mat <- SSN:::makeCovMat(dist.hydro = i.m$d, a.mat = i.m$a, b.mat = i.m$b, w.matrix = i.m$w, net.zero = n.zero, x.row = x, y.row = y, x.col = x, y.col = y, ...)
  }else{
    i.m <- getImportantMatrices.obs(d.junc, ssn@obspoints@SSNPoints[[1]]@point.data[, glmssn$args$addfunccol])
    ge <- glmssn$estimates
    ga <- glmssn$args
    cov.mat <- SSN:::makeCovMat(ge$theta, i.m$d, i.m$a, i.m$b, i.m$w, n.zero, x, y, x, y, ga$useTailDownWeight, ga$CorModels, ga$use.nugget, ga$use.anisotropy, ga$REs)
  }
  # Return
  list(
    `Covariance matrix` = cov.mat,
    `Determinant` = det(cov.mat)
  )
}