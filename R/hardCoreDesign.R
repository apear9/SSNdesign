#' The hardCoreDesign function from SSN updated to be compatible with functions from SSNDesign
#' 
#'@description
#'
#'\code{hardCoreDesign} replaces a function of the same name from the package SSN.
#' 
#'@usage
#'
#'\code{hardCoreDesign(...)} 
#'
#'@param ... Arguments for the function \code{hardCoreDesign} in the package SSN.
#'@return An object of class data.frame
#'
#'@details
#'
#'This function was written to deal with errors resulting in the \code{hardCoreDesign} function from the package SSN when it was used with SpatialStreamNetworks built from real spatial data. It is back-compatible with the \code{createSSN} function from SSN.
#' 
#' @export
hardCoreDesign <- function (n, inhibition_region, replications = 1, rep.variable = "Time", rep.values) 
{
  if (missing(rep.values)) 
    rep.values <- 1:replications
  if (replications != length(rep.values)) {
    stop("Input rep.values must contain one element for each replication")
  }
  design.function <- function(tree.graphs, edge_lengths, locations, 
                              edge_updist, distance_matrices) {
    tmp <- n
    if (length(n) == 1) 
      tmp <- rep(n, length(tree.graphs))
    if (length(inhibition_region) == 1) 
      inhibition_region <- rep(inhibition_region, length(tree.graphs))
    if (length(tree.graphs) != length(tmp) & length(tmp) != 1) {
      stop("Dimension mismatch: Input n must contain one number, or one number for each network")
    }
    if (length(tree.graphs) != length(inhibition_region) & length(inhibition_region) != 1) {
      stop("Dimension mismatch: Input inhibition_region must contain one number, or one number for each network")
    }
    initial_points <- binomialDesign(tmp)(tree.graphs, edge_lengths, locations, edge_updist, distance_matrices)
    final_result <- vector(mode = "list", length = length(initial_points))
    cumulative_locID <- 0
    for(i in 1:length(initial_points)){
      network_points.i <- initial_points[[i]]
      rids.i <- as.numeric(names(edge_updist[[i]]))
      edge_dictionary <- matrix(
        data = c(
          1:length(rids.i),
          sort(unique(rids.i))
        ),
        ncol = 2,
        byrow = FALSE
      )
      colnames(edge_dictionary) <- c("iGraph", "RealNetwork")
      ind <- match(network_points.i$edge, edge_dictionary[, 2])
      network_points.i$edge <- edge_dictionary[ind, 1]
      edgelist.i <- get.edgelist(tree.graphs[[i]])
      related_pairs <- collapseAndReorderEdgelist(edgelist.i)
      related_pairs_edges <- get.edge.ids(tree.graphs[[i]], related_pairs)
      vertex_edge_dictionary <- cbind(edgelist.i, related_pairs_edges)
      edge_lengths.i <- edge_lengths[[i]]
      edge_updist.i <- edge_updist[[i]]
      distance_matrix.i <- distance_matrices[[i]]
      colnames(distance_matrix.i) <- rownames(distance_matrix.i) <- 1:ncol(distance_matrix.i)
      network_points_edge_lengths <- edge_lengths.i[network_points.i$edge]
      rids.unique <- vertex_edge_dictionary[, 3]
      upstream_edges <- match(vertex_edge_dictionary[rids.unique,  1], vertex_edge_dictionary[,  2])
      downstream_edges <- match(vertex_edge_dictionary[rids.unique,  2], vertex_edge_dictionary[,  1])
      relationship_table_full <- cbind(vertex_edge_dictionary, upstream_edges, downstream_edges)
      colnames(relationship_table_full) <- c("UpVertex", "DownVertex", "Edge", "Is_Downstream_Of", "Is_Upstream_Of")
      dsind <- match(network_points.i$edge, relationship_table_full[, "Edge"])
      downstream_rids_per_point <- relationship_table_full[dsind, "Is_Upstream_Of"]
      point_index <- 1
      distances <- vector(mode = "numeric", nrow(network_points.i))
      while(point_index <= nrow(network_points.i)){
        current_rid <- network_points.i$edge[point_index]
        dsind <- match(current_rid, relationship_table_full[, "Edge"])
        downstream_current_rid <- relationship_table_full[dsind, "Is_Upstream_Of"]
        first <- distance_matrix.i[downstream_rids_per_point, downstream_current_rid] + network_points_edge_lengths * network_points.i$ratio + 
          network_points_edge_lengths[point_index] * network_points.i$ratio[point_index]
        second <- distance_matrix.i[downstream_rids_per_point, current_rid] + network_points_edge_lengths * network_points.i$ratio + 
          network_points_edge_lengths[point_index] * (1 - network_points.i$ratio[point_index])
        third <- distance_matrix.i[network_points.i$edge, downstream_current_rid] + network_points_edge_lengths * (1 - network_points.i$ratio) + 
          network_points_edge_lengths[point_index] * network_points.i$ratio[point_index]
        distances <- pmin(first, second, third, na.rm = TRUE)
        subset.bool <- network_points.i$edge == current_rid
        distances[subset.bool] <- network_points_edge_lengths[subset.bool] * abs(network_points.i$ratio[point_index] - network_points.i$ratio[subset.bool])
        to.keep <- distances > inhibition_region[i]
        if (any(is.na(to.keep)))
          stop("Internal error")
        to.keep[point_index] <- TRUE
        network_points.i <- network_points.i[to.keep, ]
        downstream_rids_per_point <- downstream_rids_per_point[to.keep]
        network_points_edge_lengths <- network_points_edge_lengths[to.keep]
        point_index <- point_index + 1
      }
      ind <- match(network_points.i$edge, edge_dictionary[, "iGraph"])
      network_points.i$edge <- edge_dictionary[ind, "RealNetwork"]
      network_points.i$locID <- 1:nrow(network_points.i) + cumulative_locID
      cumulative_locID <- cumulative_locID + nrow(network_points.i)
      final_result[[i]] <- network_points.i
    }
    return(final_result)
  }
  return(SSN:::replication.function(design.function, replications, 
                              rep.variable, rep.values))
}
