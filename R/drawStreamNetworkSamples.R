#' Construct probability-based and heuristic designs from Som et al. (2014)
#' 
#' @description 
#' 
#' This function is a wrapper for \code{Stream.Network.Samples}. This function constructs any of the designs outlined in Som et al. (2014), such as GRTS and stream network designs with clusters around confluences.
#' 
#' @usage 
#' 
#' \code{drawStreamNetworkSamples(ssn, new.ssn.path, overwrite.path = FALSE, sample.method, sample.size, ...)}
#' 
#' @param ssn An object of class SpatialStreamNetwork.
#' @param new.ssn.path A path for the new .ssn directory where the results should be written out.
#' @param overwrite.path A logical indicating whether the ssn directory referred to in new.ssn.path should be overwritten, if it exists already. Defaults to FALSE.
#' @param sample.method A character vector providing the label for the specific sampling design that is desired, with references to Som et al. (2014) including simple random sample "SRS", "GRTS", "GRTSmouth", "GRTSclus", H1 "Extreme.Clust.and.Singles", C3 "Trib.Sets.Ext.Singles.sample", C4 "Trib.Sets.Ext.Singles.Mouth.sample"
#' @param sample.size A numeric scalar specifying the desired sample size.
#' @param ... Other arguments to Stream.Network.Samples; chiefly, cluster.number, which is required for some sample.method arguments.
#' @return An object of class SpatialStreamNetwork. Note, any prediction points will have to be imported separately.  
#' 
#' @details 
#' 
#' STUFF.
#' 
#' @export
drawStreamNetworkSamples <- function(
  ssn, 
  new.ssn.path,
  overwrite.path = FALSE,
  sample.method,
  sample.size,
  ...
){
  
  # Check inputs
  opts <- c(
    "SRS",
    "GRTS", 
    "GRTSmouth", 
    "GRTSclus", 
    "Extreme.Clust.and.Singles", 
    "Trib.Sets.Ext.Singles.sample", 
    "Trib.Sets.Ext.Singles.Mouth.sample"
  )
  if(!is.logical(overwrite.path)){
    stop("The argument overwrite.path must be a logical.")
  }
  if(!is.character(sample.method)){
    stop("The argument sample.method must be a string.")
  }
  if(!(sample.method %in% opts)){
    msg <- paste("Invalid sampling method specified. These are your options:", paste(opts, collapse = " // "))
  }
  if(!is.numeric(sample.size)){
    stop("The argument sample.size must be a positive whole number.")
  }
  
  # Get ssn list
  ssn.list <- createSSNList(ssn)
  
  # Use Stream.Network.Samples
  result.interim <- suppressWarnings(
    Stream.Network.Samples(
      ssn.list,
      sample.method,
      sample.size,
      ...
    )
  )
  
  # Clean up output
  to.delete <- ncol(
    result.interim$ssn.obj.samples@obspoints@SSNPoints[[1]]@point.data
  )
  result.interim$ssn.obj.samples@obspoints@SSNPoints[[1]]@point.data <- result.interim$ssn.obj.samples@obspoints@SSNPoints[[1]]@point.data[,-to.delete]
  
  # Write out new ssn
  writeSSN(
    result.interim$ssn.obj.samples,
    new.ssn.path,
    o.write = overwrite.path
  )
  
  # Reimport from new path
  result.final <- importSSN(
    new.ssn.path,
    o.write = overwrite.path
  )
  
  # Return result
  return(result.final)
  
}