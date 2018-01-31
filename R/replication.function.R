replication.function <- function (design.function, replications = 1, rep.variable = "Time", 
          rep.values) 
{
  replicated.function <- function(tree.graphs, edge_lengths, 
                                  locations, edge_updist, distance_matrices) {
    n_networks <- length(tree.graphs)
    unreplicated <- design.function(tree.graphs, edge_lengths, 
                                    locations, edge_updist, distance_matrices)
    replicated <- unreplicated
    if (replications > 1) {
      for (netid in 1:n_networks) {
        for (i in 1:(replications - 1)) {
          replicated[[netid]] <- rbind(replicated[[netid]], 
                                       unreplicated[[netid]])
        }
        replicated[[netid]][, rep.variable] <- rep(rep.values, 
                                                   each = nrow(unreplicated[[netid]]))
        replicated[[netid]]$locID <- rep(unreplicated[[netid]]$locID, 
                                         times = replications)
      }
    }
    return(replicated)
  }
  return(replicated.function)
}