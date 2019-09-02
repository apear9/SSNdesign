#' Convert all components of a SpatialStreamNetwork object into \code{data.frame} objects
#' 
#' @description
#' 
#' This function takes a SpatialStreamNetwork and returns a list of one to three data.frame objects. The first data.frame object is for the stream edges and the second is for the observed points. The third data.frame is only produced when \code{preds = TRUE}. This function is intended to simplify the process of preparing a SpatialStreamNetwork to plot with ggplot2, which only accepts data.frame objects.
#' 
#' @param ssn An object of class SpatialStreamNetwork. 
#' @param preds A logical to specify whether the prediction points in the SSN (if any) should be coerced to data.frame. Defaults to FALSE.
#' @return A list. The list contains at least two elements. The first element will be the data.frame equivalent of the edges shapefile, obtained in the same format as the output of a \code{ggplot2::fortify()} call on a SpatialLines* object. This data.frame will also include the attributes of the edges. The second element will be data.frame equivalent of the obspoints slot in the SpatialStreamNetwork, including both the point.data and the point.coords slots. The third element will only be returned if \code{preds = TRUE}. This is the data.frame equivalent of the predpoints slot, processed in exactly the same way as the obspoint slot discussed previously.
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#' # Create a SSN
#' 
#' s <- createSSN(100, binomialDesign(10), 
#' binomialDesign(10), paste(tempdir(), "s.ssn", sep = "/"), TRUE)
#' 
#' # Coerce to data.frame
#' 
#' s_df <- ssn2DataFrames(ssn = s, preds = TRUE)
#' 
#' # Check elements
#' 
#' names(s_df)
#' 
#' 
#' }
#' 
#' @export
ssn2DataFrames <- function(ssn, preds = FALSE){
  
  # Input checking
  if(class(ssn) != "SpatialStreamNetwork"){
    stop("ssn must be an object of class SpatialStreamNetwork")
  }
  is.obs <- length(ssn@obspoints@SSNPoints) > 0
  if(preds){
    # Check prediction points are actually present
    is.prd <- length(ssn@predpoints@SSNPoints) > 0
    if(!is.prd){
      stop("You have set preds = TRUE but this SpatialStreamNetwork contains no prediction points.")
    }
  }
  
  if(is.obs){
    # Extract prediction and observed frames, including coordinates
    obs.cd <- ssn@obspoints@SSNPoints[[1]]@point.coords
    obs.cd <- data.frame(obs.cd)
    names(obs.cd) <- c("x", "y")
    obs.df <- getSSNdata.frame(ssn)
    obs.df <- cbind(obs.cd, obs.df)
  }
  if(preds){
    prd.cd <- ssn@predpoints@SSNPoints[[1]]@point.coords
    prd.cd <- data.frame(prd.cd)
    names(prd.cd) <- c("x", "y")
    prd.df <- getSSNdata.frame(ssn, "preds")
    prd.df <- cbind(prd.cd, prd.df)
  }
  
  # Extract data for edges
  edges.sp <- ggplot2::fortify(ssn)
  edges.df <- ssn@data
  edges.df$rid <- as.character(edges.df$rid)
  edges.df <- dplyr::left_join(edges.sp, edges.df, by = c("id" = "rid"))
  names(edges.df)[1:2] <- c("x", "y")
  
  # Assemble into list
  df.list <- list(
    edges = edges.df
  )
  if(is.obs){
    df.list$obspoints <- obs.df
  }
  if(preds){
    df.list$predpoints <- prd.df
  }
  
  # Return this list
  return(df.list)
  
}
