#' Plot the diagnostics output by the design functions
#'
#'@description 
#'
#' Creates a traceplot of the maximum utilities 
#' 
#'@usage 
#' 
#' \code{plotTrace(ssndesign)}
#' 
#'@param ssndesign The list object output by one of the main design functions.
#'@return A ggplot object
#'
#'@details 
#'
#'This function plots the trace for the utility function as the main design functions iterate through the greedy exchange algorithm. 
#'
#'Note that this function simply makes use of the output from ssndesign2DF in a standard format. Users can simply use \code{ssndesign2DF} to obtain the same data for plotting, which they can use in whatever fashion they wish.
#'  
#'@export
plotTrace <- function(ssndesign){
  
  # Check the input
  if(!is.list(ssndesign)){
    stop("The argument ssndesign must be the list output by one of the main design functions.")
  }
  
  ## Get data frame
  plot.df <- ssndesign2DF(ssndesign)
  plot.df$network <- paste("Network", plot.df$network)
  plot.df$K <- paste("K =", plot.df$K)
  
  ## Plot
  
  plot1 <- ggplot(
    data = plot.df, 
    aes(x = itr)
  ) + 
    geom_point(
      col = "red",
      aes(y = max)
    ) +
    geom_line(
      lty = 2,
      aes(y = median)
    ) +
    geom_ribbon(
      alpha = 0.5,
      fill = "grey50",
      aes(
        ymin = min,
        ymax = max
      )
    ) +
    ylab(
      "Utility value"
    ) +
    xlab(
      "Greedy exchange number"
    ) +
      theme_bw()
    
    if(length(unique(plot.df$network)) > 1){
      plot2 <- plot1 + 
        facet_grid(K ~ network)
    } else {
      plot2 <- plot1 +
        facet_wrap(~ K)
    }
  
    return(plot2)
    
}
