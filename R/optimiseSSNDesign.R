#' A function to optimise designs under pseudo-Bayesian utility functions using the Greedy Exchange Algorithm
#' 
#'@description 
#' 
#' The function \code{optimiseSSNDesign} is the main workhorse function of the package \code{SSNdesign}. It works to construct optimal and adaptive designs using a greedy exchange algorithm under a variety of different circumstances.
#' 
#'@param ssn An object of class SpatialStreamNetwork. This object must contain all potential sampling sites in its\code{obspoints} slot. 
#'@param new.ssn.path A path to a folder where the result can be stored. The folder name must end in \code{.ssn}. This function will throw an error if the argument matches the name of an existing folder. 
#'@param glmssn A fitted \code{glmssn} object. However, this model does not need to have been fitted to the \code{SpatialStreamNetwork} object in the argument ssn.
#'@param n.points A numeric or a named numeric vector specifying the size of the final design(s). See Details for more information.
#'@param legacy.sites A vector of the pids or locIDs of any 'legacy sites' which must appear in the final design. This argument is optional and should be left missing if not required. Note that the total number of sampling sites in the final design will still be \code{n.points}. 
#'@param select.by A string argument which specifies whether each pid represents a single design point or each locID. The options are "auto" (the function will autodetect whether there are multiple pids per locID and use locIDs as unique design points), "pid" or "locID". This argument defaults to "auto".
#'@param utility.function A function with the signature 'utility.function'. See Details for more information.
#'@param prior.parameters A list of functions or a matrix. If this argument is a list, the elements of this list specify independent priors for the covariance parameters in the glmssn object. If this argument is a matrix, the matrix contains the prior draws from any kind of prior (multivariate or independent) on the covariance parameters in the \code{glmssn} object. Note, in this case, the \code{n.draws} argument will be ignored because the matrix will have as many rows as there are prior draws. See Details for more information.
#'@param n.cores The number of CPUs which should be used when running \code{optimiseSSNDesign}. This argument must agree with the argument \code{parallelism}. For example, this means that, if n.cores > 1 and \code{parallelism = "none"}, this argument will be ignored and all computations will be performed sequentially. Defaults to 1. 
#'@param parallelism Must be one of "none", "windows", or "osx/linux". These arugments are insensitive to case. Note the argument must be selected appropriately for the operating system on the user's computer.  
#'@param parallelism.seed Either a numeric integer or NULL. This argument can be used to seed a random number generator which ensures reproducible calculations. This argument is effective regardless of whether parallel computations are being used.
#'@param n.optim The number of times the Greedy Exchange Algorithm is iterated to find an optimal design. Any integer greater than or equal to 1 is permissible, though larger values will produce more reliable results. Defaults to 5.
#'@param n.draws The number of Monte Carlo draws used to approximate the expected utility from the utility function per Muller (1999). Any values larger than 1 are permitted, though a minimum of 100 are recommended for the most stable results. 
#'@param extra.arguments A list of miscellaneous arguments and values which may be used to control the behaviour of the utility.function. See Details for more information.
#'@param verbose Whether messages indicating the function's progress should be printed to the console. Defaults to \code{TRUE}.
#'@param record.designs Whether the function should return all the designs evaluated. This is FALSE by default and should remain \code{FALSE} for larger examples due to the incredible amount of memory this requires.
#'@param ... Any additional arguments for the \code{foreach} iterator. The version of \code{foreach} from the package \code{doRNG} is used.
#'@return A list with the following elements: 
#' 1. ssn.old, the original and unaltered ssn. 
#' 2. ssn.new, the original ssn modified such that it contains only the sites in the optimal design.
#' 3. call, the code block used to run \code{optimiseSSNDesign}.
#' 4. glmssn, the original and unmodified glmssn.
#' 5. legacy.sites, a vector containing the locID or pid values of any legacy sites.
#' 6. utility.function, the original and unmodified utility function used to optimise the design. 
#' 7. prior.draws, a matrix containing the simulated values of the covariance parameters drawn from the priors.
#' 8. n.draws, a numeric indicating the number of Monte Carlo draws that were used when evaluating the expected utility of each design.
#' 9. random.seed, the random seed specified in the \code{parallelism.seed} argument.
#' 10. final.points, a vector of the locID or pid values for the sites included in the final design.
#' 11. best.designs.per.K.iteration, the same as final.points for each of the random starts.
#' 12. trace.per.random.start, 
#' 13. designs, a list of all the designs evaluated when optimising the design. This will be empty unless \code{record.designs = TRUE}.
#' 14. utilities, a list of all the expected utility values of the designs evaluated when optimising the design.
#' 
#'@details 
#' 
#' There are several aspects of this function that require explanations beyond the brief descriptions of the arguments above. 
#' \itemize{
#' \item n.points: this argument can be a single number or a vector. If the user supplies a single number, then, regardless of the number of isolated networks present in the ssn object, a total of n.points sites will be selected across ALL these networks. That is, the networks are not treated as separate design problems. However, if the user supplies a vector, things become more complicated. Firstly, the vector must be named. The names of the elements correspond to networks in the ssn argument and the element itself is the number of sites which should be selected within that network. The n.points argument may therefore look like this: \code{c("1" = 5, "3" = 6)}. In this example, the user is asking for 5 sites to be chosen from network 1 and 6 from network 3. It also shows that there is no need to select sites in every network. One final caveat is that no sites will appear in ssn.new for any networks which are skipped over.
#' \item utility.function: there are a number of pre-defined utility functions such as \code{DOptimality} and \code{KOptimality}. See THIS PAGE for an exhaustive list. Instructions for creating user-defined utility functions can also be found through that link.
#' \item prior.parameters: the number of elements in the prior.parameters list must be the same as the number of covariance parameters in the fitted \code{glmssn} which is also passed to this function. For example, if the only covariance parameter is the nugget (i.e. there is no spatial autocorrelation), then prior.parameters would be a list with one element. Each list element must be structured as an anonymous function as follows: \code{prior.parameters[[i]] <- function(x) runif(x)}. The random sample function does not have to be \code{runif}. The key message is that the function must only have a single argument, which is the number of draws to be taken from a random sample function. The function \code{\link{constructLogNormalPriors}} is able to construct lists of log-normal priors from \code{glmssn} objects. 
#' \item extra.arguments: this is an argument that is mostly exploited internally by \code{optimiseSSNDesign}, since the list structure is convenient for storing myriad objects such as design and distance matrices involved in the computation of the expected utility. However, some utility functions, such as \code{CPOptimality}, have additional parameters that they expect from the extra.arguments list. In particular, \code{CPOptimality} expects extra.arguments$h, a step-size argument that it uses in the forward-finite differencing of the covariance matrix. (It needs to do this to estimate the expected Fisher information matrix, which the utility function relies on.)
#' }
#' @examples
#' 
#' \dontrun{
#' ## Simulate an ssn
#' s <- createSSN(c(100, 100), binomialDesign(c(25, 25)), path = paste(tempdir(), "s1.ssn", sep = "/"), importToR = TRUE)
#' createDistMat(s)
#' ## Simulating data
#' s <- SimulateOnSSN(
#'  s, 
#'  getSSNdata.frame(s),
#'  formula = ~ 1, 
#'  coefficients = 1, 
#'  CorParms = c(1, 2, 1, 2, 1, 2, .1),
#'  addfunccol = "addfunccol"
#' )$ssn.object
#' ## Model-fitting
#' m <- glmssn(Sim_Values ~ 1, s, addfunccol = "addfunccol")
#' ## Construct a list of log-normal priors
#' p <- constructLogNormalPriors(m1)
#' ## Use optimiseSSNDesign
#' r.together <- optimiseSSNDesign(
#'  ssn = s, new.ssn.path = paste(tempdir(), "s2.ssn", sep = "/"), glmssn = m, n.points = 25, 
#'  utility.function = DOptimality, prior.parameters = p, n.cores = 2, parallelism = "windows", 
#' parallelism.seed = 123, n.optim = 1, n.draws = 100)
#' r.apart <- optimiseSSNDesign(
#'  ssn = s, new.ssn.path = paste(tempdir(), "s3.ssn", sep = "/"), glmssn = m, n.points = c("1" = 7, "2" = 8), 
#'  utility.function = DOptimality, prior.parameters = p, n.cores = 2, parallelism = "windows", 
#'  parallelism.seed = 123, n.optim = 1, n.draws = 100)
#' }
#' 
#'@references
#'
#' Muller, P. (1999). Simulation Based Optimal Design. \emph{Bayesian Statistics}, \emph{6}, 1-13.   
#' 
#'@export
optimiseSSNDesign <- function(
  ssn, 
  new.ssn.path,
  glmssn, 
  n.points, 
  legacy.sites,
  select.by = "auto",
  utility.function,
  prior.parameters,
  n.cores = 1,
  parallelism = "none",
  parallelism.seed = NULL,
  n.optim = 5, 
  n.draws = 500,
  extra.arguments = list(),
  verbose = TRUE,
  record.designs = FALSE,
  ...
){
  
  # Early exits for bad inputs
  if(!isSSN(ssn)){
    stop("The argument SSN must be of class SpatialStreamNetwork.")
  }
  if(dir.exists(new.ssn.path)){
    stop("The SSN folder specified in the argument new.ssn.path already exists. Please choose a new one or delete the existing folder.")
  }
  if(!is.character(parallelism)){
    stop("The parallelism must be one of the following strings, with any capitalisation: none, windows, osx/linux.")
  }
  parallelism <- stringr::str_to_lower(parallelism)
  if(!(parallelism %in% c("none", "windows", "osx/linux"))){
    stop("The argument parallelism must be one of : none, windows, osx/linux.")
  }
  n.draws <- floor(n.draws)
  if(n.draws < 1){
    stop("The argument n.draws must be a positive integer greater than 1. Please make sure this is the case. Note that the recommended minimum value for n.draws is 100.")
  }
  if(!(is.null(parallelism.seed) | is.numeric(parallelism.seed))){
    stop("The argument parallelism.seed should either be left as NULL or changed to a positive integer.")
  }
  if(is.numeric(parallelism.seed)){
    if(parallelism.seed <= 0){
      stop("The argument parallelism.seed must be a POSITIVE integer.")
    }
    if(floor(parallelism.seed) != parallelism.seed){
      stop("The argument parallelism.seed must be an INTEGER such that floor(parallelism.seed) is equal to parallelism.seed")
    }
  }
  
  # Extract a bunch of information we will eventually return
  the.call <- match.call()
  
  # Check if any points are to remain fixed
  any.fixed  <- !missing(legacy.sites)
  
  # Detect parallel options
  if(n.cores > 1 & parallelism == "none"){
    stop("If the argument n.cores > 1, then the argument parallelism must be either ")
  }
  
  # Force parallelism = "none" if only one core to be used
  if(n.cores == 1){
    parallelism <- "none" # no point in setting up a cluster, esp. if it can cause issues depending on OS.
  }
  
  # Set cluster type to PSOCK if windows OS, FORK otherwise. 
  # If no parallelism req'd, then set to run sequentially.
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
  
  # Set parallel seed if present
  if(!is.null(parallelism.seed)){
    registerDoRNG(parallelism.seed)
  } # Works the same way for registerDoSEQ
  
  # Extract out ssn data frame
  ssn.obs.data <- getSSNdata.frame(ssn)
  
  # Detect whether we're dealing with PIDs or locIDs
  if(select.by == "auto"){
    by.locID <- length(unique(ssn.obs.data$locID)) != length(ssn.obs.data$pid)
  } else if(select.by == "locID"){
    by.locID <- TRUE
    if(verbose) message("The argument select.by has been set to 'locID'. The resulting design will be for points with unique locIDs. Please note that this means the arguments n.points and legacy.sites will be interpreted as referring to pids also.")
  } else {
    by.locID <- FALSE
  }
  
  # If locIDs, map locIDs to PIDs
  if(by.locID){
    if(verbose & select.by == "auto"){
      message("Multiple pids per locID have been detected. The resulting design will be for points with unique locIDs.")
    }
    alls <- ssn.obs.data$locID
    lIDs <- unique(alls)
    pIDs <- ssn.obs.data$pid
    locIDs <- vector("list", length(lIDs))
    for(i in 1:length(lIDs)){
      locIDs[[i]] <- pIDs[alls == lIDs[i]]
    }
    names(locIDs) <- lIDs
  }
  
  # Find number of points in each network
  netIDs <- ssn.obs.data$netID
  nets <- unique(netIDs)
  n.nets <- length(nets)
  p.nets <- vector("numeric", n.nets)
  cnt <- 1
  for(i in nets){
    p.nets[cnt] <- sum(netIDs == i)
    cnt <- cnt + 1
  }
  names(p.nets) <- nets
  
  # Extract out max number of points
  d.nets.names <- names(n.points)
  if(length(n.points) != 1){
    ind <- tryCatch(match(d.nets.names, names(p.nets)), error = function(e) stop("n.points must be a named vector, with the names of the elements corresponding to netIDs which exist in the SSN object."))
    stopifnot(all(unlist(lapply(1:length(n.points), function(x, ind){n.points[x] < p.nets[ind[x]]}, ind = ind))))
  }
  
  # Subtract the number of fixed points from the total number of points requested
  if(any.fixed){
    if(length(n.points) == 1){
      n.points <- n.points - length(legacy.sites)
      if(n.points < 1){
        stop("The number of legacy sites exceeds the number of sites requested for the final design.")
      }
    } else {
      # This is a bit painful if there are multiple networks...
      # 1. Find the networks for which designs have been requested
      # 2. Find the number of fixed points belonging to each network
      # 3. Subtract and check that all are greater than 0
      for(i in 1:length(n.points)){
        d.net <- names(n.points)[i]
        ssn.data.i <- ssn.obs.data[which(ssn.obs.data$netID == d.net), ]
        if(by.locID){
          n.f.in.net <- sum(legacy.sites %in% ssn.data.i$locID)
        } else {
          n.f.in.net <- sum(legacy.sites %in% ssn.data.i$pid)
        }
        n.points[i] <- n.points[i] - n.f.in.net
      }
      if(any(n.points < 1)){
        stop("At least one network has too many legacy sites.")
      }
    }
  }
  
  # Extract out distance matrices and stuff
  d <- getStreamDistMatInOrder(ssn)
  d <- constructTotalMatrix(d)
  n.z <- d$net.zero
  afv <- glmssn$args$addfunccol # additive function value for constructing matrices
  if(is.null(afv)){
    afv.col.o <- rep(1, nrow(ssn@obspoints@SSNPoints[[1]]@point.data))
  } else {
    afv.col.o <- ssn.obs.data[, afv]
  }
  M <- getImportantMatrices.obs(d$d.junc, afv.col.o) # M for Many matrices hehe
  
  # Put these matrices in the extra.arguments list in the way expected by the utility functions
  Matrices.Obs <- list(Matrices.Obs = M, net.zero.obs = n.z)
  extra.arguments <- append(
    extra.arguments, 
    Matrices.Obs
  )
  
  # Extract out model matrix using the observed data and the formula from the glmssn object
  m.form <- glmssn$args$formula
  m.char <- as.character(m.form)
  m.form <- as.formula(paste(m.char[1], m.char[3]))
  # if(m.char[3] != "1"){
  #   obs.X  <- SSN:::dataXY(m.form, ssn.obs.data, glmssn$args$family, glmssn$args$trialscol, glmssn$args$trans.power, glmssn$args$trans.shift, order(ssn.obs.data$pid), glmssn$args$CorModels)$Xmats$X2
  # } else {
  #   obs.X <- model.matrix(m.form, data = ssn.obs.data)
  # }
  obs.X <- model.matrix(m.form, data = ssn.obs.data)
  obs.C  <- ssn@obspoints@SSNPoints[[1]]@point.coords
  row.names(obs.X) <- row.names(obs.C) <- ssn.obs.data$pid
  colnames(obs.C) <- c("x", "y")
  extra.arguments$obs.X <- obs.X
  extra.arguments$obs.C <- obs.C
  
  # It may be important to do the same thing for the prediction points, so we detect whether they exist and if so we construct the matrices.
  if(anyPreds(ssn)){
    # Two sets of matrices -- one for pred x pred and another for pred x obs
    d <- getStreamDistMatInOrder(ssn, "preds")
    d <- constructTotalMatrix(d)
    n.z <- d$net.zero
    if(is.null(afv)){
      afv.col.p <- rep(1, nrow(ssn@predpoints@SSNPoints[[1]]@point.data))
    } else {
      afv.col.p <- getSSNdata.frame(ssn, "preds")[, afv]
    }
    M <- getImportantMatrices.obs(d$d.junc, afv.col.p)
    # Add these to the extra.arguments list
    Matrices.prd <- list(Matrices.prd = M, net.zero.prd = n.z)
    # The second set of matrices requires a slightly different set of code
    d <- getStreamDistMatInOrder.predsxobs(ssn)
    d <- constructTotalMatrix(d, pxo = TRUE)
    n.z <- d$net.zero
    M <- getImportantMatrices.pxo(d$d.junc,t(d$d.junc),afv.col.o,afv.col.p)
    # Add to extra.arguments list
    Matrices.pxo <- list(Matrices.pxo = M, net.zero.pxo = n.z)
    # Extraction of model matrix
    prd.X <- model.matrix(m.form, getSSNdata.frame(ssn, "preds"))
    prd.C <- ssn@predpoints@SSNPoints[[1]]@point.coords
    row.names(prd.X) <- row.names(prd.C) <- getSSNdata.frame(ssn, "preds")$pid
    colnames(prd.C) <- c("x", "y")
    # Now stick these into the extra.arguments list
    extra.arguments <- append(extra.arguments, Matrices.prd)
    extra.arguments <- append(extra.arguments, Matrices.pxo)
    extra.arguments$prd.X <- prd.X
    extra.arguments$prd.C <- prd.C
  }
  
  # Take prior draws now, and fix these for later
  prior.draws <- matrix()
  if(length(prior.parameters) != 0){
    if(!is.matrix(prior.parameters)){
      n.parms <- length(prior.parameters)
      prior.draws <- matrix(0, ncol = n.parms, nrow = n.draws)
      for(parameter in 1:n.parms){
        prior.draws[, parameter] <- prior.parameters[[parameter]](n.draws)
      }
    } else {
      prior.draws <- prior.parameters
      n.draws <- nrow(prior.draws)
    }
    
  }
  # As required for ED and EK optimality
  fep <- MASS::mvrnorm(n.draws, glmssn$estimates$betahat, glmssn$estimates$covb) # The fixed effects
  fep <- unname(fep)
  extra.arguments$Empirical.FEP <- fep
  
  # Construct sets of candidate points
  if(by.locID){
    if(length(n.points) != 1){
      c.points <- vector("list", length(n.points))
      for(i in 1:length(n.points)){
        ind <- ssn.obs.data$netID == d.nets.names[i]
        c.points[[i]] <- ssn.obs.data$locID[ind]
      }
    } else {
      c.points <- list(ssn.obs.data$locID)
    }
  } else {
    if(length(n.points) != 1){
      c.points <- vector("list", length(n.points))
      for(i in 1:length(n.points)){
        ind <- ssn.obs.data$netID == d.nets.names[i]
        c.points[[i]] <- ssn.obs.data$pid[ind]
      }
    } else {
      c.points <- list(ssn.obs.data$pid)
    }
  }
  
  # Convert all of these to numeric
  for(q in 1:length(c.points)){
    c.points[[q]] <- anum(c.points[[q]])
  }
  
  # If legacy.sites not missing, do the same for the fixed sites. Remove them from the pool of points which may be selected.
  if(any.fixed){
    x.points <- vector("list", length(n.points))
    # Two cases:
    # 1. All networks are part of a single large design problem
    # 2. Each network is its own design problem
    if(length(n.points) == 1){
      # Just remove all the fixed points from the c.points list
      ind <- !(c.points[[1]] %in% legacy.sites)
      x.points[[1]] <- legacy.sites
      c.points[[1]] <- c.points[[1]][ind]
    } else {
      # For each network which is being designed for, remove the corresponding points from legacy.sites
      for(i in 1:length(n.points)){
        ind1 <- !(c.points[[i]] %in% legacy.sites)
        ind2 <- !ind1
        x.points[[i]] <- c.points[[i]][ind2]
        c.points[[i]] <- c.points[[i]][ind1]
      }
    }
    # Then if each fixed point represents a locID, then convert these to PIDs.
    if(by.locID){
      x.points.eval <- vector("list", length(n.points))
    } else {
      x.points.eval <- x.points
    }
  } else {
    # Generate blank list elements to append
    x.points <- x.points.eval <- lapply(1:length(n.points), function(x) numeric(0))
  }
  
  # Loop through to optimise
  f.p.nets <- vector("list", length(n.points))
  f.u.nets <- vector("list", length(n.points))
  for(i in 1:length(n.points)){
    
    n.final <- n.points[i]
    c.point <- c.points[[i]]
    x.point <- x.points[[i]]
    x.point.eval <- x.points.eval[[i]]
    best.u.in.opt <- c()
    best.in.optim <- list()
    us.per.optims <- list()
    designs.per.optims <- list()
    trace.per.optims <- list()
    for(j in 1:n.optim){
      
      # Record time
      begin.time <- Sys.time()
      
      # Print progress to user
      if(verbose){
        message(paste0("Now on iteration ", j, " out of ", n.optim))
      }
      
      # Take random sample of points
      r.point <- sample(c.point, n.final, FALSE)
      if(by.locID){
        r.point.eval <- c()
        for(z in 1:n.final){
          ind <- match(r.point[z], lIDs)
          r.point.eval <- c(r.point.eval, locIDs[[ind]])
        }
        r.point.eval <- c(r.point.eval, x.point.eval)
      } else {
        r.point.eval <- c(r.point, x.point.eval) 
      }
      
      # Evaluate utility once, store information
      u <- u.stars <- u.star <- max.utility <- utility.function(ssn, glmssn, r.point.eval, prior.draws, n.draws, extra.arguments)
      u.star.old <- u.star - 0.0001 # Just to make the while loop condition true in the first instance
      r <- r.point
      
      # Initialise loop
      designs <- list(r.point)
      d.star <- designs[[1]]
      u.all <- c(u)
      while(u.star > u.star.old){
        
        u.star.old <- u.star
        
        for(k in 1:n.final){
          
          # Construct list of designs to evaluate
          u.this.run <- c()
          designs.this.run <- designs.this.run.eval <- list()
          
          # Identify candidate points to exchange
          ind <- !(c.point %in% r.point) 
          not.in <- c.point[ind]
          
          # Replace k'th position with by cycling through remaining values
          for(l in 1:length(not.in)){
            # Coordinate exchange
            r.point[k] <- not.in[l]
            # locID-pid matching if required
            if(by.locID){
              r.point.eval <- c()
              for(z in 1:n.final){
                ind <- match(r.point[z], lIDs)
                r.point.eval <- c(r.point.eval, locIDs[[ind]])
              }
              r.point.eval <- c(r.point.eval, x.point.eval)
            } else {
              r.point.eval <- c(r.point, x.point.eval)
            }
            # Put these in a list of designs
            designs.this.run <- append(designs.this.run, list(r.point))
            designs.this.run.eval <- append(designs.this.run.eval, list(r.point.eval))
          }
          
          # Evaluate expected utilities here
          n.eval <- length(designs.this.run.eval)
          u.this.run <- foreach(x = 1:n.eval, .combine = c, .packages = c("SSN", "SSNdesign"), ...) %dopar% {
            utility.function(ssn, glmssn, designs.this.run.eval[[x]], prior.draws, n.draws, extra.arguments)
          }
          designs <- append(designs, designs.this.run)
          u.all <- c(u.all, u.this.run)
          max.utility <- max(u.this.run)
          # Check logical condition
          if(max.utility > u.star){
            better.ind <- which(u.this.run == max.utility)[1]
            # ind <- which(u.this.run == u.star)
            d.star <- designs.this.run[[better.ind]]
            r.point <- d.star
            u.star <- max.utility
          }
          u.stars <- c(u.stars, u.star)
        }
      }
      
      # Store best
      best.u.in.opt <- c(best.u.in.opt, u.star)
      best.in.optim <- append(best.in.optim, list(c(d.star, x.point)))
      us.per.optims <- append(us.per.optims, list(u.all))
      
      # Store designs that were evaluated 
      if(record.designs) designs.per.optims <- append(designs.per.optims, list(designs))
      
      # Save algorithm trace
      trace.per.optims <- append(trace.per.optims, list(u.stars))
      
      # Record end time and print to user
      end.time <- Sys.time()
      if(verbose){
        elapsed.time <- end.time - begin.time
        message(paste0("Iteration ", j, " took ", round(as.numeric(elapsed.time), 2), " ", units(elapsed.time)))
      }
      
    }
    
    ind <- which(best.u.in.opt == max(best.u.in.opt))
    f.p.nets[[i]] <- best.in.optim[[ind[1]]]
    f.u.nets[[i]] <- us.per.optims
    
  }
  
  # Collapse results into one list
  f.p <<- unlist(f.p.nets)
  f.p.scope <- f.p
  
  # Subset based on list
  if(by.locID){
    ssn.new <- subsetSSN(ssn, new.ssn.path, locID %in% f.p)
  } else {
    ssn.new <- subsetSSN(ssn, new.ssn.path, pid %in% f.p)
  }
  
  # Delete the object that was assigned in the global environment
  rm(list = "f.p", pos = ".GlobalEnv")
  
  # Create return list
  if(!record.designs) designs.per.optims <- NULL
  if(missing(legacy.sites)) legacy.sites <- numeric(0)
  ssn.design <- list(
    ssn.old = ssn, 
    ssn.new = ssn.new, 
    call = the.call,
    glmssn = glmssn,
    legacy.sites = legacy.sites,
    utility.function = utility.function,
    prior.draws = prior.draws,
    n.draws = n.draws,
    random.seed = parallelism.seed,
    final.points = f.p.scope, 
    best.designs.per.K.iteration = best.in.optim, 
    trace.per.random.start = trace.per.optims,
    designs = designs.per.optims,
    utilities = f.u.nets
  )
  class(ssn.design) <- "ssndesign"
  
  # Return list, subset, and diagnostics
  return(
    ssn.design
  )
  
}