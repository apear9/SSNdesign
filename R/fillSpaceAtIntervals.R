fillSpaceAtIntervals <- function(d.m, spacing, row = 1){
  
  # d.m is a distance matrix
  # row is the row of the distance matrix that remains fixed for the space-filling design
  if(row < 1){
    stop("Row must be a positive integer.")
  }
  if(floor(row) != ceiling(row)){
    stop("Row must be a whole number")
  }
  if(spacing <= 0){
    stop("spacing must be a positive number")
  }
  
  the.row <- d.m[row, ]
  
  if(all(d.m == 0)){
    return(colnames(d.m))
  } # In case all sites at the same place
  
  # Set up spacings
  # Find max distance
  max.dist <- max(the.row)
  spaces <- seq(0, max.dist, spacing)
  nspaces <- length(spaces)
  
  # For each 'space', find an appropriate point or point(s) at that distance from the 0 element in the distance matrix
  selected <- c() # vector to store pids
  for(i in 1:nspaces){
    selected <- c(selected, closestToZero(spaces[i], the.row))
  }
  
  # Return vector of pids
  return(unique(selected))
  
}
