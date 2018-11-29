#' Evaluate fixed designs
#' 
#'@description 
#'
#' Evaluates the utility of a list of user-specified designs. 
#' 
#'@usage 
#' 
#' \code{evaluateFixedDesigns(ssn, glmssn, afv.column, list.designs, utility.function, prior.parameters, n.draws = 500, extra.arguments = NULL)}
#' 
#'@param ssn an object of class SpatialStreamNetwork
#'@param glmssn an object of class glmssn
#'@param afv.column the name of the column in the SpatialStreamNetwork object that contains the additive function values
#'@param list.designs a list containing vectors of  
#'@param utility.function a function with the signature utility.function. Users may define their own. This package provides several built-in utility functions. 
#'@param prior.parameters a function to act as a prior for covariance parameter values
#'@param n.draws a numeric value for the number of Monte Carlo draws to take when evaluating potential designs
#'@param extra.arguments a list of any extra parameters which can be used to control the behaviour of this function or the utility function
#'@return A dataframe of design ID (integers from 1 in the same order as list.designs); Size, the number of design points; and, Utility, the utility value. 
#'
#'@details 
#'
#'This function takes a list of user-specified designs for a spatial stream network and evaluates some utility function over those designs. This function should be relatively fast compared to the functions which search for an optimal design. However, because this function is not parallellised, it is recommended that users limit themselves to fewer than one hundred designs. Note, this function is called internally in the design diagnostic functions.  
#'  
#'@export
evaluateFixedDesigns <- function(
  ssn,
  glmssn, 
  afv.column, 
  list.designs, 
  utility.function, 
  prior.parameters, 
  n.draws = 500, 
  extra.arguments = NULL
){
  # Check inputs
  
  if(!is(ssn, "SpatialStreamNetwork")){
    stop("The argument 'ssn' must be an object of class SpatialStreamNetwork.")
  }
  n <- nnetwork(ssn)
  if(!is.numeric(n.draws) | n.draws != floor(n.draws)){
    stop("The argument n.draws must be a whole-numbered numeric value.")
  }
  if(!is.numeric(n.draws) | n.draws < 1){
    stop("The argument n.draws must have a positive integer value.")
  }
  if(n.draws < 2){
    stop("Please choose a sensible number of draws for n.draws. Recommended values are 100 or more.")
  }
  
  ## Find number of designs
  nd <- length(list.designs)
  
  ## Find the length of each design
  ndp <- lapply(list.designs, length)
  size.offsets <- rep(0, nd)
  
  ### Code to detect whether there are replicates on sites
  n.pids <- length(
    unique(
      ssn@obspoints@SSNPoints[[1]]@point.data$pid
    )
  )
  n.locIDs <- length(
    unique(
      ssn@obspoints@SSNPoints[[1]]@point.data$locID
    )
  )
  
  ## Extract K for greedy exchange algorithm
  if(is.null(extra.arguments)){
    extra.arguments <- list()
  }
  if(is.null(extra.arguments$K)){
    K <- 20 
  } else {
    K <- extra.arguments$K
  }
  
  ## Extract important matrices for each network
  dist.junc.obs <- getStreamDistMatInOrder(ssn)
  if(length(ssn@predpoints@SSNPoints) > 0){
    dist.junc.prd <- getStreamDistMatInOrder(ssn, "preds")
    dist.junc.pxo <- getStreamDistMatInOrder.predsxobs(ssn)
  }
  ## Construct important matrices and shove them in to the extra.arguments list
  total.mats.obs <- constructTotalMatrix(dist.junc.obs)
  extra.arguments$net.zero.obs <- total.mats.obs$net.zero
  # Process matrix to obtain important matrices
  extra.arguments$Matrices.Obs <- getImportantMatrices.obs(
    total.mats.obs$d.junc,
    ssn@obspoints@SSNPoints[[1]]@point.data[, afv.column]
  )
  rm(total.mats.obs) # remove from memory in case it is very large
  # Do the same for prediction related matrices if needed
  if(length(ssn@predpoints@SSNPoints) > 0){
    total.mats.prd <- constructTotalMatrix(dist.junc.prd)
    extra.arguments$net.zero.prd <- total.mats.prd$net.zero
    extra.arguments$Matrices.prd <- getImportantMatrices.obs(
      total.mats.prd$d.junc,
      ssn@predpoints@SSNPoints[[1]]@point.data[, afv.column]
    )
    rm(total.mats.prd)
    total.mats.pxo <- constructTotalMatrix(dist.junc.pxo, TRUE)
    extra.arguments$net.zero.pxo <- total.mats.pxo$net.zero
    extra.arguments$Matrices.pxo <- getImportantMatrices.pxo(
      total.mats.pxo$d.junc,
      t(total.mats.pxo$d.junc),
      ssn@obspoints@SSNPoints[[1]]@point.data[, afv.column],
      ssn@predpoints@SSNPoints[[1]]@point.data[, afv.column]
    )
    rm(total.mats.pxo)
  }
  
  # Initialise loops by creating a list of potential sites to choose from
  is.replicated <- n.pids != n.locIDs
  if(is.replicated){
    print("Replicates found. Mapping PIDs to locIDs...")
    # Create mapping from locID to pid
    all.locIDs <- as.numeric(as.character(ssn@obspoints@SSNPoints[[1]]@point.data$locID))
    unique.locIDs <- unique(all.locIDs)
    locID.to.pid <- vector("list", n.locIDs)
    # Get all the PIDs associated with each locID
    for(i in 1:n.locIDs){
      locID.i <- all.locIDs == i
      locID.to.pid[[i]] <- ssn@obspoints@SSNPoints[[1]]@point.data$pid[locID.i]
    }
    # Map all designs expressed in terms of locIDs to PIDs
    for(i in 1:nd){
      ndp.i <- ndp[[i]]
      design.to.eval <- c()
      for(j in 1:ndp.i){
        ind <- which(unique.locIDs == list.designs[[i]][j])
        if(length(ind) < 1){
          size.offsets[i] <- size.offsets[i] + 1
          next
        }
        design.to.eval <- c(design.to.eval, locID.to.pid[[ind]])
      }
      list.designs[[i]] <- design.to.eval
    }
  }
  
  ## Set up dataframe to store results
  m <- matrix(
    nrow = nd, 
    ncol = 3  # ID, Size, and Utility 
  ) 
  m[, 1] <- 1:nd 
  m[, 2] <- unlist(ndp) - size.offsets # sometimes sites in the designs aren't actually present among the locIDs...
  
  ## Evaluate all designs
  for(i in 1:nd){
    U <- utility.function(
      ssn,
      glmssn,
      list.designs[[i]],
      prior.parameters, 
      n.draws,
      extra.arguments
    )
    m[i, 3] <- U
  }
  
  # Turn into real data frame and name columns
  results.df <- data.frame(m)
  names(results.df) <- c("ID", "Size", "Utility") 
  
  # Possibly rename elements in first column
  if(!is.null(names(list.designs))){
    results.df[, 1] <- names(list.designs)
  }
  
  # Spit out the results
  return(results.df)
  
}
