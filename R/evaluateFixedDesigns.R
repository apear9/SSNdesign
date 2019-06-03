#' Evaluate fixed designs
#' 
#'@description 
#'
#' Evaluates the utility of a list of user-specified designs. 
#' 
#'@param ssn an object of class SpatialStreamNetwork
#'@param glmssn an object of class glmssn
#'@param list.designs a list containing vectors of pid or locID values. Use pid values if the sites in a design change over time.
#'@param list.of a string (either "pid" or "locID") indicating whether the designs are expressed as vectors of pid or locID values
#'@param utility.function a function with the signature utility.function. Users may define their own. This package provides several built-in utility functions. 
#'@param prior.parameters a function to act as a prior for covariance parameter values
#'@param n.draws a numeric value for the number of Monte Carlo draws to take when evaluating potential designs. This defaults to 1000.
#'@param extra.arguments a list of any extra parameters which can be used to control the behaviour of this function or the utility function
#'@return A \code{data.frame} with the columns: ID (integers from 1 in the same order as list.designs); Size, the number of design points; Expected Utility; Efficiency, the ratio between each expected utility and the largest expected utility; and Efficiency_Unlogged, the same but for expected utilities expressed on the log scale. 
#'
#'@details 
#'
#' This function takes a list of user-specified designs for a spatial stream network and evaluates some utility function over those designs. This function should be relatively fast compared to the functions which search for an optimal design. However, because this function is not parallellised, it is recommended that users limit themselves to fewer than one hundred designs. Note, this function is called internally in the design diagnostic functions.  
#'  
#'@export
evaluateFixedDesigns <- function(
  ssn,
  glmssn, 
  list.designs, 
  list.of = "pid",
  utility.function, 
  prior.parameters, 
  n.draws = 1000, 
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
  
  ## Get additive function values
  afv.column <- glmssn$args$addfunccol
  
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
  
  ## Extract model frame and coordinates of the observed points
  ssn.obs.data <- getSSNdata.frame(ssn)
  m.form <- glmssn$args$formula
  m.char <- as.character(m.form)
  m.form <- as.formula(paste(m.char[1], m.char[3]))
  obs.X  <- model.matrix(m.form, ssn.obs.data, contrasts)
  obs.C  <- ssn@obspoints@SSNPoints[[1]]@point.coords
  row.names(obs.X) <- row.names(obs.C) <- ssn.obs.data$pid
  colnames(obs.C) <- c("x", "y")
  extra.arguments$obs.X <- obs.X
  extra.arguments$obs.C <- obs.C
  
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
  if(!is.null(afv.column)){
    extra.arguments$Matrices.Obs <- getImportantMatrices.obs(
      total.mats.obs$d.junc,
      ssn@obspoints@SSNPoints[[1]]@point.data[, afv.column]
    )
  } else {
    extra.arguments$Matrices.Obs <- getImportantMatrices.obs(
      total.mats.obs$d.junc,
      rep(1, nrow(ssn.obs.data))
    )
  }
  rm(total.mats.obs) # remove from memory in case it is very large
  # Do the same for prediction related matrices if needed
  if(anyPreds(ssn) > 0){
    total.mats.prd <- constructTotalMatrix(dist.junc.prd)
    extra.arguments$net.zero.prd <- total.mats.prd$net.zero
    if(!is.null(afv.column)){
      extra.arguments$Matrices.prd <- getImportantMatrices.obs(
        total.mats.prd$d.junc,
        ssn@predpoints@SSNPoints[[1]]@point.data[, afv.column]
      )
    } else {
      extra.arguments$Matrices.prd <- getImportantMatrices.obs(
        total.mats.prd$d.junc,
        rep(1, nrow(getSSNdata.frame(ssn, "preds")))
      )
    }
    rm(total.mats.prd)
    total.mats.pxo <- constructTotalMatrix(dist.junc.pxo, TRUE)
    extra.arguments$net.zero.pxo <- total.mats.pxo$net.zero
    if(!is.null(afv.column)){
      extra.arguments$Matrices.pxo <- getImportantMatrices.pxo(
        total.mats.pxo$d.junc,
        t(total.mats.pxo$d.junc),
        ssn@obspoints@SSNPoints[[1]]@point.data[, afv.column],
        ssn@predpoints@SSNPoints[[1]]@point.data[, afv.column]
      )
    } else {
      extra.arguments$Matrices.pxo <- getImportantMatrices.pxo(
        total.mats.pxo$d.junc,
        t(total.mats.pxo$d.junc),
        rep(1, nrow(ssn.obs.data)),
        rep(1, nrow(getSSNdata.frame(ssn, "preds")))
      )
    }
    rm(total.mats.pxo)
    # Similarly extract the required model matrices and coordinates for the prediction points
    prd.X <- model.matrix(m.form, getSSNdata.frame(ssn, "preds"))
    prd.C <- ssn@predpoints@SSNPoints[[1]]@point.coords
    row.names(prd.X) <- row.names(prd.C) <- getSSNdata.frame(ssn, "preds")$pid
    colnames(prd.C) <- c("x", "y")
    # Now stick these into the extra.arguments list
    extra.arguments$prd.X <- prd.X
    extra.arguments$prd.C <- prd.C
  }
  
  # Take prior draws now, and fix these for later
  prior.draws <- matrix()
  if(length(prior.parameters) != 0){
    n.parms <- length(prior.parameters)
    prior.draws <- matrix(0, ncol = n.parms, nrow = n.draws)
    for(parameter in 1:n.parms){
      prior.draws[, parameter] <- prior.parameters[[parameter]](n.draws)
    }
  }
  # As required for ED and EK optimality
  fep <- MASS::mvrnorm(n.draws, glmssn$estimates$betahat, glmssn$estimates$covb) # The fixed effects
  fep <- unname(fep)
  extra.arguments$Empirical.FEP <- fep
  
  # Initialise loops by creating a list of potential sites to choose from
  # is.replicated <- n.pids != n.locIDs
  is.replicated <- list.of != "pid"
  if(is.replicated){
    message("Mapping PIDs to locIDs...")
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
      prior.draws, 
      n.draws,
      extra.arguments
    )
    m[i, 3] <- U
  }
  
  # Turn into real data frame and name columns
  results.df <- data.frame(m)
  names(results.df) <- c("ID", "Size", "Expected utility") 
  
  # Possibly rename elements in first column
  if(!is.null(names(list.designs))){
    results.df[, 1] <- names(list.designs)
  }
  
  # Compute efficiencies
  max.exp <- max(results.df$`Expected utility`)
  results.df$Efficiency <- results.df$`Expected utility`/max.exp
  results.df$Efficiency_Unlogged <- exp(results.df$`Expected utility` - max.exp)
  
  # Spit out the results
  return(results.df)
  
}
