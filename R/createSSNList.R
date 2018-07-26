createSSNList <- function(ssn){
  
  # Get networks with sites on them
  nets <- unique(ssn@obspoints@SSNPoints[[1]]@point.data$netID)
  
  # Return ssn.list
  list(
    ssn.object = ssn,
    bin.table = getBIDtables(ssn, nets)
  )
  
}
