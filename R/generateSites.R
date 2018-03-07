#' Simulate observed and/or predicted points on an empty SpatialStreamNetwork.
#' 
#' \code{generateSites()} simulates observed and predicted sites on an empty SpatialStreamNetwork
#' 
#' @param ssn an object of class SpatialStreamNetwork
#' @param edgeweights a string; the name of a column in the data slot of the SpatialStreamNetwork that is used to weight edges
#' @param obsDesign a design function
#' @param predDesign a design function 
#' @param importToR a logical value
#' @param treeFunction DEPRECATED
#' @return An object of class SpatialStreamNetwork
#' 
#' @examples 
#' \dontrun{NONE AS YET}
#' 
#' @section Warning:
#' This function is currently written such that any observed or predicted sites generated on a SpatialStreamNetwork will overwrite any existing sites.
#' 
#' @export
generateSites <- function (
  ssn, 
  edgeweights,
  edgeafv,
  obsDesign, 
  predDesign = noPoints,
  importToR = FALSE, 
  treeFunction = igraphKamadaKawai
){
  path = ssn@path
  if (missing(obsDesign)) 
    stop("Input obsDesign cannot be missing")
  if (missing(path)) 
    stop("Path cannot be missing")
  if (missing(ssn)) {
    stop("Input ssn cannot be missing")
  }
  if (length(path) != 1) 
    stop("Please enter a single path")
  info <- file.info(path)
  isdir <- info$isdir
  if (is.na(isdir)) {
    dir.create(path)
  }
  else if (isdir == FALSE) {
    stop("Unable to create directory")
  }
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(path)
  if(is.character(ssn)){
    ssn <- importStreams(ssn)
  }
  if(is.null(edgeweights)){
    edgeweights <- "shreve"
    edgeafv <- "addfunccol"
  }
  if(is.null(edgeafv) & edgeweights != "shreve"){
    stop("Please specify an existing column to serve as the set of additive function values on this SSN. Alternatively, leave edgeweights = NULL to generate additive function values based on Shreve's stream order.")
  }
  n_networks <- nnetwork(ssn)
  # Reformat netId and networkID columns to be numeric
  ssn@data$netID <- as.numeric(as.character(ssn@data$netID))
  ssn@network.line.coords$NetworkID <- as.numeric(as.character(ssn@network.line.coords$NetworkID))
  edges <- vector(mode = "list", length = n_networks)
  tree.graphs <- edges
  locations <- edges
  rids <- edges
  edge_lengths <- edges
  edge_updist <- edges
  line_data <- data.frame()
  # empty dist matrix, list to contain n_network matrices
  drvr <- dbDriver("SQLite")
  conn <- dbConnect(drvr, paste(ssn@path, "binaryID.db", sep = "/"))
  distance_matrices <- list()
  # for loop over each network
  for (netid in 1:n_networks) {
    #subnetwork <- subset(ssn, netID == netid)
    edges_this_network <- subset(ssn@network.line.coords, NetworkID == netid)
    edge_updist_this_network <- edges_this_network$DistanceUpstream
    #edge_updist_this_network <- edge_updist[[netid]] 
    #edges_this_network <- edges[[netid]] # extract edges belonging to current network in loop
    #binary_id_table <- binary_ids_tables[[netid]] # extract table of netIDs belonging to current network in loop
    binary_id_table <- dbReadTable(conn, paste0("net", netid))
    binary_id_table <- binary_id_table[order(as.numeric(binary_id_table$rid)), ]
    distance_matrix <- matrix(0, nrow(binary_id_table), nrow(binary_id_table)) # zero matrix on network
    colnames(distance_matrix) <- rownames(distance_matrix) <- binary_id_table$rid
    partial_match_function <- function(binary_id1, binary_id2) {
      min_len <- min(nchar(binary_id1), nchar(binary_id2))
      for (j in 1:min_len) {
        if (substr(binary_id1, j, j) != substr(binary_id2,j, j)) 
          return(j - 1)
      }
      return(min_len)
    }
    # Deriving network distance by partial string matching
    character_binary_ids <- as.character(binary_id_table$binaryID)
    for (i in 1:nrow(binary_id_table)) {
      current_binary_id <- binary_id_table$binaryID[i]
      current_rid <- as.character(binary_id_table$rid[i])
      current_updist <- edges_this_network$DistanceUpstream[as.character(edges_this_network$SegmentID) == current_rid]
      matching_characters <- sapply(character_binary_ids, 
                                    partial_match_function, as.character(current_binary_id))
      matching_substring <- substr(binary_id_table$binaryID, 1, matching_characters)
      indices <- match(matching_substring, binary_id_table$binaryID)
      if (any(is.na(indices))) 
        stop("Internal Error")
      downstream_rids <- binary_id_table$rid[indices]
      network_rids <- rownames(edges_this_network)
      indices <- matchIndices(downstream_rids, network_rids)
      downstream_updists <- edges_this_network$DistanceUpstream[indices]
      distance_matrix[as.character(current_rid), ] <- pmax(current_updist - downstream_updists, rep(0, length(downstream_updists)))
    }
    reindex <- match(binary_id_table$rid, edges_this_network$SegmentID)
    distance_matrix <- distance_matrix[reindex, reindex]
    distance_matrices[[netid]] <- distance_matrix + t(distance_matrix)
  }
  dbDisconnect(conn)
  # derive shreve stream orders if no edgeweights and edgeafv values are provided
  if(edgeweights == "shreve"){
    conn <- dbConnect(drvr, paste(ssn@path, "binaryID.db", sep = "/"))
    for(i in 1:n_networks){
      bin.id.tab <- dbReadTable(conn, paste0("net", i))
      if(i == 1){
        shreve.frame <- calculateShreveStreamOrder(bin.id.tab)
        mshreve <- max(shreve.frame$shreve)
        shreve.frame[, edgeafv] <- shreve.frame$shreve/mshreve
      } else {
        temp.frame <- calculateShreveStreamOrder(bin.id.tab)
        mshreve <- max(temp.frame$shreve)
        temp.frame[, edgeafv] <- temp.frame$shreve/mshreve
        shreve.frame <- rbind(shreve.frame, temp.frame)
      }
    }
    delete.old.shreve <- grep(edgeweights, names(ssn@data))
    delete.old.afv <- grep(edgeafv, names(ssn@data))
    delete.old <- c(delete.old.shreve, delete.old.afv)
    if(length(delete.old) > 0){
      ssn@data <- ssn@data[, -delete.old] 
    }
    ssn@data <- merge(ssn@data, shreve.frame, by = "rid")
    dbDisconnect(conn)
  }
  graphinfo <- readshpnw(ssn)
  graphs <- nel2igraph(graphinfo[[2]], graphinfo[[3]], eadf = graphinfo[[5]], Directed = TRUE)
  line_data <- data.frame(
    rid = edge_attr(graphs)["rid"],
    netID = edge_attr(graphs)["netID"],
    weights = edge_attr(graphs)[edgeweights],
    addfunccol = edge_attr(graphs)[edgeafv]
  )
  names(line_data)[3:4] <- c(edgeweights, edgeafv)
  for(i in 1:n_networks){
   # print(i)
    tree.graphs[[i]] <- subgraph.edges(graphs, which(E(graphs)$netID == i))
    attributes.i <- E(tree.graphs[[i]])
    edge_lengths[[i]] <- attributes.i$Shape_Leng
    rids[[i]] <- attributes.i$rid
    edge_updist[[i]] <- attributes.i$upDist
    names(edge_updist[[i]]) <- rids[[i]]
    locations[[i]] <- layout.auto(tree.graphs[[i]])
    edges[[i]] <- get.edgelist(tree.graphs[[i]])
  #  print(rids[[i]])
  }
  ## Using the designs to generate observed and predicted sites
  obs_sites <- obsDesign(tree.graphs, edge_lengths, locations,
                         edge_updist, distance_matrices) 
  pred_sites <- predDesign(tree.graphs, edge_lengths, locations,
                           edge_updist, distance_matrices)
 # print(obs_sites)
  # character.order <- as.character(order(as.character(1:n_networks)))
  # numeric.order <- order(as.numeric(character.order))
  max_observed_locID <- max(unlist(lapply(obs_sites, function(x) max(x$locID))))
  for (i in 1:length(pred_sites)) {
    pred_sites[[i]]$locID <- pred_sites[[i]]$locID + max_observed_locID
  }
  n_obs_sites <- unlist(lapply(obs_sites, function(x) return(dim(x)[1])))
  n_pred_sites <- unlist(lapply(pred_sites, function(x) return(dim(x)[1])))
  sites_data <- data.frame()
  combined_site_location_data <- c()
  pred_data <- data.frame()
  combined_pred_location_data <- c()
  cumulative_pids <- 0
  ## Something goes on here: locating all points
  f <- function(row){
    rid <- as.integer(row[1])
    ratio <- as.numeric(row[2])
    location <- locatePointOnEdge(ssn, rid, ratio) # break point???
    ret <- c(
      location[1:2], 
      edge_updist_this_network[names(edge_updist_this_network) == as.character(rid)] - location[3]
    )
    return(ret)
  }
  for (netid in 1:n_networks) {
    edge_lengths_this_network <- edge_lengths[[netid]]
    edges_this_network <- edges[[netid]]
    edge_updist_this_network <- edge_updist[[netid]]
    locations_this_network <- locations[[netid]]
    n_locations_this_network <- n_obs_sites[netid] + n_pred_sites[netid]
    rids_this_network <- rids[[netid]]
    pred_sites_this_network <- pred_sites[[netid]]
    obs_sites_this_network <- obs_sites[[netid]] 
    obs_location_data <- data.frame(
      rid = obs_sites_this_network$edge,
      ratio = obs_sites_this_network$ratio, 
      locID = obs_sites_this_network$locID,
      stringsAsFactors = FALSE
    )
    # print(edge_updist_this_network)
    # print(obs_sites_this_network) ## break point?
    # print(obs_location_data)      ## break point?
    pred_location_data <- data.frame(
      rid = pred_sites_this_network$edge,
      ratio = pred_sites_this_network$ratio,
      locID = pred_sites_this_network$locID,
      stringsAsFactors = FALSE
    )
    # print("HERE")
    if (n_locations_this_network > 0) {
      pred_location_data_this_network <- t(apply(pred_location_data, 1, f))
      if (length(pred_location_data_this_network) > 0) {
        colnames(pred_location_data_this_network) <- c("NEAR_X","NEAR_Y","upDist")
      }
      else {
        pred_location_data_this_network <- matrix(0, 0, 3)
        colnames(pred_location_data_this_network) <- c("NEAR_X","NEAR_Y","upDist")
      }
      obs_location_data_this_network <- t(apply(obs_location_data, 1, f))
      if (length(obs_location_data_this_network) > 0) {
        colnames(obs_location_data_this_network) <- c("NEAR_X","NEAR_Y","upDist")
      } else {
        obs_location_data_this_network <- matrix(0, 0, 3)
        colnames(obs_location_data_this_network) <- c("NEAR_X","NEAR_Y","upDist")
      }
      if (n_obs_sites[netid] > 0) {
        obs_pids <- (1:n_obs_sites[netid]) + cumulative_pids
      }
      else obs_pids <- integer(0)
      if (n_pred_sites[netid] > 0) {
        pred_pids <- n_obs_sites[netid] + (1:n_pred_sites[netid]) +
          cumulative_pids
      }
      else pred_pids <- integer(0)
    }
    else {
      obs_location_data_this_network <- pred_location_data_this_network <- data.frame(
        NEAR_X = numeric(0),
        NEAR_Y = numeric(0), 
        upDist = numeric(0)
      )
      obs_pids <- integer(0)
      pred_pids <- integer(0)
    }
    cumulative_pids <- cumulative_pids + n_locations_this_network
    obs_data_this_network <- data.frame(
      locID = obs_location_data[, "locID"], 
      upDist = obs_location_data_this_network[, "upDist"], 
      pid = obs_pids,
      netID = rep(netid, length(obs_pids)),
      rid = obs_location_data[, "rid"], 
      ratio = obs_location_data[, "ratio"], 
      weights = line_data[match(obs_location_data[, "rid"], line_data[, "rid"]), edgeweights], 
      addfunccol = line_data[match(obs_location_data[, "rid"], line_data[, "rid"]), edgeafv], 
      stringsAsFactors = FALSE
    )
    if (ncol(obs_sites_this_network) > 3) {
      obs_data_this_network <- cbind(
        obs_data_this_network,
        obs_sites_this_network[, -match(c("edge", "ratio", "locID"), colnames(obs_sites_this_network)), drop = FALSE]
      )
    }
    rownames(obs_data_this_network) <- obs_pids
    rownames(obs_location_data_this_network) <- obs_pids
    pred_data_this_network <- data.frame(
      locID = pred_location_data[, "locID"],
      upDist = pred_location_data_this_network[, "upDist"],
      pid = pred_pids,
      netID = rep(netid, length(pred_pids)),
      rid = pred_location_data[, "rid"], 
      ratio = pred_location_data[, "ratio"], 
      weights = line_data[match(pred_location_data[, "rid"], line_data[, "rid"]), edgeweights], 
      addfunccol = line_data[match(pred_location_data[, "rid"], line_data[, "rid"]), edgeafv], stringsAsFactors = FALSE)
    if (ncol(pred_sites_this_network) > 3) {
      pred_data_this_network <- cbind(
        pred_data_this_network,
        pred_sites_this_network[, -match(c("edge", "ratio", "locID"), colnames(pred_sites_this_network)), drop = FALSE]
      )
    }
    rownames(pred_data_this_network) <- pred_pids
    rownames(pred_location_data_this_network) <- pred_pids
    if (n_obs_sites[netid] > 0) {
      sites_data <- rbind(obs_data_this_network[1:n_obs_sites[netid],
                                                , drop = FALSE], sites_data)
      combined_site_location_data <- rbind(
        obs_location_data_this_network[1:n_obs_sites[netid], , drop = FALSE], 
        combined_site_location_data
      )
    }
    if (n_pred_sites[netid] > 0) {
      pred_data <- rbind(pred_data_this_network[(1:n_pred_sites[netid]), , drop = FALSE], pred_data)
      combined_pred_location_data <- rbind(
        pred_location_data_this_network[(1:n_pred_sites[netid]), , drop = FALSE],
        combined_pred_location_data
      )
    }
  }
  # Write sites and predictions to file
  if (length(combined_site_location_data) == 0){
    stop("At least one observation site must be present")
  }
  sites <- SpatialPointsDataFrame(combined_site_location_data[, c("NEAR_X", "NEAR_Y"), drop = FALSE], sites_data, match.ID = TRUE)
  writeOGR(sites, ".", "sites", verbose = FALSE, driver = "ESRI Shapefile", overwrite_layer = TRUE)
  if (length(combined_pred_location_data) > 0) {
    preds <- SpatialPointsDataFrame(combined_pred_location_data[,
                                                                c("NEAR_X", "NEAR_Y"), drop = FALSE], pred_data,
                                    match.ID = TRUE)
    writeOGR(preds, ".", "preds", verbose = FALSE, driver = "ESRI Shapefile", overwrite_layer = TRUE)
  }
  ind1 <- colnames(sites@data) == c("netID")
  ind2 <- colnames(sites@data) == c("rid")
  ind3 <- colnames(sites@data) == c("upDist")
  if (sum(ind1) == 0) {
    stop("netID is missing from sites attribute table")
  }
  if (sum(ind2) == 0) {
    stop("rid is missing from sites attribute table")
  }
  if (sum(ind3) == 0) {
    stop("upDist is missing from sites attribute table")
  }
  if (is.factor(sites@data$netID)) {
    sites@data$netID <- as.character(sites@data$netID)
  }
  network.point.coords <- data.frame(sites@data[, "netID"], 
                                     sites@data[, "rid"], sites@data[, "upDist"])
  colnames(network.point.coords) <- c("NetworkID", "SegmentID", "DistanceUpstream")
  network.point.coords <- as.data.frame(network.point.coords)
  row.names(network.point.coords) <- row.names(sites@data)
  attributes(network.point.coords)$locID <- as.numeric(as.character(sites@data$locID))[sites@data$locID]
  network.point.coords[, 1] <- as.factor(network.point.coords[, 1])
  network.point.coords[, 2] <- as.factor(network.point.coords[, 2])
  network.point.coords[, 3] <- as.numeric(network.point.coords[, 3])
  rm(ind1, ind2, ind3)
  op <- new("SSNPoint", network.point.coords = network.point.coords, 
            point.coords = sites@coords, point.data = sites@data, 
            points.bbox = sites@bbox)
  ssn@obspoints@SSNPoints[[1]] <- op
  ssn@obspoints@ID[[1]] <- "Obs"
  rm(network.point.coords, sites, op)
  if (length(combined_pred_location_data) > 0) {
    rownames(preds@data) <- preds@data[, "pid"]
    rownames(preds@coords) <- preds@data[, "pid"]
    preds@data$locID <- as.factor(preds@data$locID)
    if (is.factor(preds@data$netID)) {
      preds@data$netID <- as.character(preds@data$netID)
    }
    network.point.coords <- data.frame(preds@data[, "netID"], preds@data[, "rid"], preds@data[, "upDist"])
    colnames(network.point.coords) <- c("NetworkID", "SegmentID", "DistanceUpstream")
    network.point.coords <- as.data.frame(network.point.coords)
    row.names(network.point.coords) <- row.names(preds@data)
    attributes(network.point.coords)$locID <- as.numeric(levels(preds@data$locID))[preds@data$locID]
    network.point.coords[, 1] <- as.factor(network.point.coords[, 1])
    network.point.coords[, 2] <- as.factor(network.point.coords[, 2])
    pp <- new("SSNPoint", network.point.coords = network.point.coords, 
              point.coords = preds@coords, point.data = preds@data, 
              points.bbox = preds@bbox)
    ssn@predpoints@SSNPoints[[1]] <- pp
    ssn@predpoints@ID[[1]] <- "preds"
    rm(preds, pp, network.point.coords)
  }
  setwd(old_wd)
  if (importToR) {
    if (sum(n_pred_sites) > 0)
      return(importSSN(path, predpts = "preds", o.write = TRUE))
    else return(importSSN(path, o.write = TRUE))
  }
  else return(ssn)
}
