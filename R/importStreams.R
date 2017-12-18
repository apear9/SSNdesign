#' Import a shapefile of stream edges without observed or predicted sites as a SpatialStreamNetwork
#' 
#' \code{importStreams()} imports a SpatialStreamNetwork that has no observed or predicted sites.
#' 
#' @param filepath a path to a .ssn folder
#' @param o.write a logical value
#' @return An object of class SpatialStreamNetwork
#' 
#' @examples
#' \dontrun{NONE AS YET}
#' 
#' @export
importStreams <- function(filepath, o.write = FALSE){
  old.wd <- getwd()
  Path <- dirname(filepath)
  ssn.obj <- basename(filepath)
  if (Path == ".") {
    Path = old.wd
  }
  setwd(paste(Path, "/", ssn.obj, sep = ""))
  edges <- readOGR(".", "edges", verbose = FALSE, stringsAsFactors = FALSE, 
                   integer64 = "allow.loss")
  rownames(edges@data) <- edges@data[, "rid"]
  if (exists("edges") == 0) {
    stop("edges.shp is missing from ", Path, " folder")
  }
  if (getinfo.shape("edges.shp")[[2]] != 3 & getinfo.shape("edges.shp")[[2]] != 
      13) {
    stop("edges.shp does not have polyline geometry")
  }
  ind1 <- colnames(edges@data) == c("netID")
  ind2 <- colnames(edges@data) == c("rid")
  ind3 <- colnames(edges@data) == c("upDist")
  if (sum(ind1) == 0) {
    stop("netID is missing from streams attribute table")
  }
  if (sum(ind2) == 0) {
    stop("rid is missing from streams attribute table")
  }
  if (sum(ind3) == 0) {
    stop("upDist is missing from streams attribute table")
  }
  if (is.factor(edges@data$netID)) {
    edges@data$netID <- as.character(edges@data$netID)
  }
  network.line.coords <- data.frame(edges@data$netID, edges@data[, 
                                                                 "rid"], edges@data[, "upDist"])
  colnames(network.line.coords) <- c("NetworkID", "SegmentID", 
                                     "DistanceUpstream")
  network.line.coords <- as.data.frame(network.line.coords)
  row.names(network.line.coords) <- row.names(edges@data)
  network.line.coords[, 1] <- as.factor(network.line.coords[, 
                                                            1])
  network.line.coords[, 2] <- as.factor(network.line.coords[, 
                                                            2])
  rm(ind1, ind2, ind3)
  op <- new("SSNPoint", network.point.coords = data.frame(NULL), 
            point.coords = matrix(nrow = 1, ncol=2), point.data = data.frame(NULL), 
            points.bbox = matrix(nrow = 2, ncol = 2))
  ops <- new("SSNPoints")
  ops@SSNPoints[[1]] <- op
  ops@ID[[1]] <- "Obs"
  ssn <- new("SpatialStreamNetwork", edges, network.line.coords = network.line.coords, 
             obspoints = ops, 
             #predpoints = pps, 
             path = paste(Path, 
                                                             "/", ssn.obj, sep = ""))
  ssn@obspoints@SSNPoints[[1]]@point.data$netID <- as.factor(ssn@obspoints@SSNPoints[[1]]@point.data$netID)
  ssn@data$netID <- as.factor(ssn@data$netID)
  rm(network.line.coords, edges)
  createBinaryID(ssn, o.write = o.write)
  setwd(old.wd)
  return(ssn)
}
