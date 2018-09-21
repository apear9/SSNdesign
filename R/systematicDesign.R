#' The systematicDesign function from SSN updated to be compatible with functions from SSNDesign
#' 
#'@description
#'
#'\code{systematicDesign} replaces a function of the same name from the package SSN.
#' 
#'@usage
#'
#'\code{systematicDesign(...)} 
#'
#'@param ... Arguments for the function \code{systematicDesign} as in the package SSN.
#'@return An object of class data.frame.
#'
#'@details
#'
#'This function was written to deal with errors resulting in the \code{systematicDesign} function from the package SSN when it was used with SpatialStreamNetworks built from real spatial data. It is back-compatible with the \code{createSSN} function from SSN.
#' 
#' @export
systematicDesign <- function(spacing, replications = 1, rep.variable = "Time", rep.values) {
  if (missing(rep.values)) 
    rep.values <- 1:replications
  if (replications != length(rep.values)) {
    stop("Input rep.values must contain one element for each replication")
  }
  design.function <- function(tree.graphs, edge_lengths, locations, 
                              edge_updist, distance_matrices) {
    if (length(spacing) == 1) 
      spacing <- rep(spacing, length(tree.graphs))
    if (length(spacing) != length(tree.graphs)) {
      stop("Dimension mismatch: Input spacing must contain one number, or one number for each network")
    }
    n_networks <- length(tree.graphs)
    result <- vector(mode = "list", length = length(n_networks))
    cumulative_locID <- 0
    for (netid in 1:n_networks) {
      spacing_this_network <- spacing[netid]
      graph <- tree.graphs[[netid]]
      edge_lengths_this_network <- edge_lengths[[netid]]
      rids <- names(edge_updist[[netid]])
      edges_this_network <- get.edgelist(graph)
      points_this_network <- sort(unique(as.numeric(edges_this_network)))
      positions_per_segment <- vector(mode = "list", length = length(edge_lengths_this_network))
      done_points <- !(points_this_network %in% edges_this_network[, 2])
      done_segments <- c()
      segment_remaining <- c()
      while (length(done_segments) != nrow(edges_this_network)) {
        can_calculate <- done_points[match(edges_this_network[, 1], points_this_network)]
        can_calculate[done_segments] <- FALSE
        can_calculate_indices <- which(can_calculate)
        if (!any(can_calculate)) 
          stop("Internal error")
        for (index in can_calculate_indices) {
          edge <- edges_this_network[index, ]
          remaining <- segment_remaining[match(match(edge[1], edges_this_network[, 2]), done_segments)]
          if (is.null(remaining)) 
            remaining <- spacing_this_network
          if(is.na(remaining))
            remaining <- spacing_this_network
          edge_length <- edge_lengths_this_network[index]
          if (edge_length + remaining < spacing_this_network) {
            segment_remaining <- c(segment_remaining, 
                                   edge_length + remaining)
          }
          else {
            positions_per_segment[[index]] <- seq(spacing_this_network - 
                                                    remaining, edge_length, by = spacing_this_network)
            segment_remaining <- c(segment_remaining, 
                                   edge_length - max(positions_per_segment[[index]]))
          }
          done_segments <- c(done_segments, index)
          done_points[match(edge[2], points_this_network)] <- TRUE
        }
      }
      proportions_per_segment <- positions_per_segment
      for (i in 1:length(proportions_per_segment)) proportions_per_segment[[i]] <- proportions_per_segment[[i]]/edge_lengths_this_network[i]
      unreplicated <- data.frame(edge = rep(rids, times = unlist(lapply(proportions_per_segment, 
                                                                        length))), ratio = unlist(proportions_per_segment), 
                                 stringsAsFactors = FALSE)
      unreplicated$locID <- 1:nrow(unreplicated) + cumulative_locID
      cumulative_locID <- cumulative_locID + nrow(unreplicated)
      result[[netid]] <- unreplicated
    }
    return(result)
  }
  return(SSN:::replication.function(design.function, replications, 
                              rep.variable, rep.values))
}
