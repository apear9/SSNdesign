#' A function to calculate the Shreve stream order of a stream edge
#' 
#'@description
#'
#'For a given binaryID table, this function quickly computes the Shreve stream order of every stream segment by counting the number of upstream inlet edges there are.
#'
#'@usage
#'
#'\code{calculateShreveStreamOrder(binary.id.table)}
#'
#'@param binary.id.table A data.frame representing a binaryID table extracted from the binaryID.db in the path of a SpatialStreamNetwork.
#'@return A data.frame with two columns: the first is the RID of each stream segment; the second is the Shreve order of each segment.
#'
#'@details
#'
#'This function uses string matching on the binaryIDs devised by ... (2010) to calculate the Shreve stream order of every stream segment in a network. This process is made computationally efficient by the realisation that the Shreve stream order is the number of inlet segments upstream of a given segment.
#'
#'@export
calculateShreveStreamOrder <- function(binary.id.table){
  # Find all inlets; inlets are all edges whose binary ids do not appear in the binary ids of any edge 'above' them
  bids <- binary.id.table$binaryID
  bids <- vapply(
    bids,
    function(x){gsub(" ", "", x)},
    vector("character", 1)
  )
  bids <- unname(bids)
  ind.inlet <- vapply(
    bids,
    function(x){
      viable.bids <- bids[nchar(bids) <= nchar(x) + 1]
      match.positions.last <- stringr::str_locate(viable.bids, x)[, 2]
      viable.positions <- match.positions.last <= nchar(x)
      sum(str_count(viable.bids[viable.positions], x), na.rm = T) == 1
    },
    vector("logical", 1)
  )
  ind.inlet <- unname(ind.inlet)
  # find upstream relationships for all edges and count number of inlets upstream (i.e. shreve order)
  shreve.order <- vapply(
    bids,
    function(x){
      match.positions.last <- stringr::str_locate(bids, x)[, 2]
      viable.positions <- match.positions.last <= nchar(x)
      inlets.upstream <- viable.positions & ind.inlet
      return(sum(inlets.upstream, na.rm = TRUE))
    },
    vector("numeric", 1)
  )
  # return dataframe associating shreve orders with rids
  shreve.order.with.rids <- data.frame(
    rid = binary.id.table$rid, 
    shreve = unname(shreve.order)
  )
  names(shreve.order.with.rids) <- c("rid", "shreve")
  return(shreve.order.with.rids)
}

