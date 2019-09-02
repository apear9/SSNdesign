#' Plot trace of the Greedy Exchange Algorithm
#' 
#' @description 
#' 
#' A generic plot function adapted for `ssndesign` objects. It plots the traces from the optimisation algorithm.
#' 
#' @method plot ssndesign
#' 
#' @param x An object of class \code{ssndesign}
#' @param y Not used 
#' @param which.iteration A numeric indicating the trace that should be plotted. This can be a vector, in which case all traces corresponding to those indices will be plotted. All traces are plotted by default. 
#' @param legend Whether the legend should be drawn on the plot. Defaults to \code{FALSE}. Not advised if there are many traces because the legend will be too large to properly fit on the plot. 
#' @param ... Additional arguments to the \code{plot} function.
#' 
#' @details 
#' 
#' The plot method produces line plots tracing the expected utility of the best design in the algorithm for each set of coordinate exchanges.
#' 
#' @method plot ssndesign
#' @export plot.ssndesign
plot.ssndesign <- function(x, y, which.iteration = 1:length(x$trace.per.random.start), legend = TRUE, ...){
  # Check that no one has input the y argument
  if(!missing(y)){
    stop("Not valid for plot.ssndesign...")
  }
  
  # Check that which.iteration is numeric and has length > 0
  if(!is.numeric(which.iteration) | !(length(which.iteration) > 0)){
    stop("The argument which.iteration must be numeric.")
  }
  
  # Check that which.iteration does not exceed reasonable bounds
  if(max(which.iteration) > length(x$trace.per.random.start)){
    stop("No element of the argument which.iteration should exceed the length of the number of traces in the argument x.")
  }
  
  # Plot
  plot.data <- x$trace.per.random.start
  plot.data.max <- unlist(lapply(plot.data, max))
  plot.data.min <- unlist(lapply(plot.data, min))
  plot.data.len <- unlist(lapply(plot.data, length))
  plot.data <- plot.data[order(plot.data.max, decreasing = TRUE)]
  plot(plot.data[[which.iteration[1]]], type = "l", col = 1, xlab = paste("Exchange number (recycles every ", (length(x$final.points) - length(x$legacy.sites)), "iterations)"), ylab = "Best U(d)", ylim = c(min(plot.data.min), max(plot.data.max)), xlim = c(1, max(plot.data.len)), ...)
  if(length(which.iteration) > 1){
    for(i in 2:length(which.iteration)){
      lines(1:length(plot.data[[which.iteration[i]]]), plot.data[[which.iteration[i]]], col = i, ...)
    }
    if(legend) legend("bottomright", legend = paste("Random start", which.iteration), lty = 1, col = 1:length(which.iteration))
  }
  
}