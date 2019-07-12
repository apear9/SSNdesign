  #' Extract a likelihood profile for the covariance parameters in a \code{glmssn}
#' 
#' @description 
#' 
#' This function computes the likelihood surface in the vicinity of the maximum likelihood estimates of the covariance parameters. Note: this is experimental. This function was initially developed as a workaround to the problem of defining log-normal priors for the covariance parameters when the standard errors extracted by \code{\link{getExpectedSE}} were manifestly unreasonable. 
#' 
#' @param glmssn An object of class \code{glmssn}
#' @param which.vary A vector indicating which of the covariance parameters should be varied. You must vary at least one parameter but you can vary all of them if it suits you.
#' @param grid.parameters An OPTIONAL argument. You can pass a list giving the grid of parameter values for which to evaluate the log-likelihood surface. This can be left missing. In this case, the function attempts to guess the grid range and resolution based on the log-scale confidence-interval on each of the varying parameters.
#' @param parallelism An OPTIONAL argument. This argument indicates what operating system a parallel cluster should be set up for. This option can be one of "none" (for no parallelism), "osx/linux" or "windows". This is not case sensitive. The input should match the user's operating system. Note that this argument is ignored if the next argument \code{n.cores = 1}. 
#' @param n.cores An OPTIONAL argument. This argument specifies how many CPUs should be assigned to compute the log-likeihoods in parallel. This argument can be any positive integer.
#' @param ... Any additional arguments to \code{foreach}. Users are NOT recommended to fiddle with this. 
#' @return A \code{data.frame} with one column for each covariance parameter and another column for the log-likelihood. 
#' 
#' @author Jay Ver Hoef and Alan R. Pearse
#' 
#' @export
getLikelihoodProfile <- function(glmssn, which.vary, grid.parameters, parallelism = "none", n.cores = 1, ...){
  
  # Checks
  if(class(glmssn) != "glmssn"){
    stop("The argument glmssn must be of class glmssn.")
  }
  n.vary <- length(which.vary)
  n.parm <- length(glmssn$estimates$theta)
  if(n.vary > n.parm | n.vary == 0){
    stop("Wrong number of parameters to vary.")
  } 
  if(any(which.vary > n.parm)){
    stop("At least one element of which.vary does not reference a valid covariance parameter.")
  }
  if(n.cores > 1 & parallelism == "none"){
    stop("If the argument n.cores > 1, then the argument parallelism must be either ")
  }
  
  ## Set up parallel clusters
  parallelism <- stringr::str_to_lower(parallelism)
  
  # Force parallelism = "none" if only one core to be used
  if(n.cores == 1){
    parallelism <- "none" # no point in setting up a cluster, esp. if it can cause issues depending on OS.
  }
  
  # Set cluster type to PSOCK if windows OS, FORK otherwise. If no parallelism req'd, then set to run sequentially.
  if(parallelism == "windows"){
    cl <- makeCluster(n.cores, type = "PSOCK")
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # Kill cluster in case of crash or happy-path-exit
  } else {
    if(parallelism == "osx/linux"){
      cl <- makeCluster(n.cores, type = "FORK")
      registerDoParallel(cl)
      on.exit(stopCluster(cl)) # Kill cluster in case of crash or happy-path-exit
    }else{
      registerDoSEQ() # for plebs
    }
  }
  
  # Get an SSN object; need to evaluate profile likelihood
  ssn <- glmssn$ssn.object
  
  # Generate arguments for m2LLstream
  data_ <- ssn@obspoints@SSNPoints[[1]]@point.data
  data_ <- cbind(data_, ssn@obspoints@SSNPoints[[1]]@point.coords)
  distord <- order(data_[,"netID"],data_[,"pid"])
  names(distord) <- rownames(data_)[distord]
  xcol <- "coords.x1"
  ycol <- "coords.x2"
  family <- tolower(glmssn$args$family)
  
  # Get data in design matrix, etc.
  dataXY.out <- SSN:::dataXY(glmssn$args$formula, data_,
                             family = family, trialscol = NULL,
                             trans.power = glmssn$args$trams.power,
                             trans.shift = glmssn$args$trans.shift,
                             CorModels = glmssn$args$CorModels,
                             distord = distord)
  
  REs <- dataXY.out$REs
  REmodelmatrices <- dataXY.out$REmodelmatrices
  n.all <- dataXY.out$sampsizes$n.all
  ind <- dataXY.out$indvecs$ind.allxy
  
  xcoord <- data_[distord,xcol]
  ycoord <- data_[distord,ycol]
  xcoord.data <- xcoord[ind]
  ycoord.data <- ycoord[ind]
  
  z <- zt <- dataXY.out$respvecs$z
  X2 <- dataXY.out$Xmats$X2
  n.allxy <- dataXY.out$sampsizes$n.allxy
  
  # Create the special distance matrices
  mats <- getStreamDistMatInOrder(ssn)
  total.mats <- constructTotalMatrix(mats)
  net.zero <- total.mats$net.zero
  d.junc <- total.mats$d.junc
  i.mats <- getImportantMatrices.obs(d.junc, data_[, glmssn$args$addfunccol])
  
  # Extract out the individual matrices we need to pass as arguments to ProfLike
  dm <- i.mats$d
  am <- i.mats$a
  bm <- i.mats$b
  cm <- i.mats$c
  wm <- i.mats$w
  
  # Get estimates
  theta <- glmssn$estimates$theta
  log.theta <- log(theta[1:length(theta)])
  
  if(missing(grid.parameters)){
    # Try to guess a good grid resolution and range
    log.se <- getExpectedSE(glmssn, TRUE)
    log.lower <- log.theta - log.se * 2 # VERY APPROXIMATE
    log.upper <- log.theta + log.se * 2 # Yup
    grid.parameters <- vector("list", n.parm)
    for(i in 1:n.vary){
      i.ind <- which.vary[i]
      grid.parameters[[i.ind]] <- c(log.lower[i.ind], log.upper[i.ind])
    }
    null.ind <- unlist(lapply(grid.parameters, is.null))
    # Slot in remaining parameters as fixed values
    grid.parameters[null.ind] <- log.theta[null.ind]
  } 
  
  # Construct true grid as data.frame
  grid.cols <- vector("character", n.parm)
  fixd.cols <- unlist(lapply(grid.parameters, function(x) length(x) == 1))
  for(i in 1:n.parm){ # problem here... what if parameter DOESN'T VARY
    if(fixd.cols[i]){
      to.add <- paste0("V", i, " = ", grid.parameters[[i]])
    } else {
      to.add <- paste0("V", i, " = seq(", grid.parameters[[i]][1], ", ", grid.parameters[[i]][2],", length.out = 201)")
    }
    grid.cols[i] <- to.add
  }
  grid.code <- paste0(
    "expand.grid(", paste(grid.cols, collapse = ", "), ")"
  )
  grid.vals <- eval(parse(text = grid.code))
  
  # Now evaluate over this grid. 
  n.combs <- nrow(grid.vals)
  result.vec <- foreach(i = 1:n.combs, .combine = c, .packages = c("SSN", "SSNdesign"), ...) %dopar% {
    # log.theta.i <- unlist(grid.vals[i, ])
    ProfLike(
      unlist(grid.vals[i, ]), 
      zt = zt,
      X = X2,
      dist.hydro = dm,
      weight = wm,
      net.zero = net.zero,
      a.mat = am,
      b.mat = bm,
      x.dat = xcoord.data,
      y.dat = ycoord.data,
      CorModels = glmssn$args$CorModels,
      use.nugget = glmssn$args$use.nugget,
      use.anisotropy = glmssn$args$use.anisotropy,
      EstMeth = glmssn$args$EstMeth,
      REs = REs, 
      scale = attributes(theta)$scale,
      maxrang = NULL
    )
  }
  
  # Join loglikelihoods to the data.frame
  grid.vals$loglik <- result.vec
  
  # Return result
  grid.vals
  
}
