#' A function to optimise designs under pseudo-Bayesian utility functions using the Greedy Exchange Algorithm
#' 
#'@description 
#' 
#' \code{optimiseSSNDesign} is the main workhorse function of the package \code{SSNdesign}. It works with a range of utilities to construct optimal solutions to both static and adaptive design problems.
#' 
#'@usage
#' 
#' \code{optimiseSSNDesign(ssn, new.ssn.path, glmssn, n.points, legacy.sites, utility.function, prior.parameters, n.cores = 1, parallelism = "none", parallelism.seed = NULL, n.optim = 5, n.draws = 500, ...)}
#' 
#'@param ssn An object of class SpatialStreamNetwork. Its obspoints slot must contain all the potential sampling sites for the stream network. It must also contain all the covariates and other information referred to in the glmssn object. 
#'@param new.ssn.path A path to a folder where the result can be stored. 
#'@param glmssn An object of class glmssn. This must be a fitted model. However, this model does not need to have been fitted to the SpatialStreamNetwork data in the argument ssn.
#'@param n.points A numeric or a named numeric vector specifying the size of the final design(s). See Details for more information.
#'@param legacy.sites A numeric, character, or factor vector indicating any 'legacy sites' which must appear in the final design. This argument is OPTIONAL. Simply leave blank if not required. Note that the total number of sampling sites in the final design will still be n.points. 
#'@param utility.function A function with the signature 'utility.function'. See Details for more information.
#'@param prior.parameters A list of functions. The elements of this list specify the independent priors on the covariance parameters in the glmssn object. See Details for more information.
#'@param n.cores The number of CPUs which should be used when running \code{optimiseSSNDesign}. This argument must agree with the argument parallelism. For example, if n.cores > 1 and \code{parallelism = "none"}, this argument will be ignored and all computations will be performed sequentially. Defaults to 1. 
#'@param parallelism A character which must be one of "none", "windows", "linux/osx". These can be spelled with in any case or in any combination of cases. Note the argument must be selected appropriately for the operating system on the user's computer. Definitely do not select "linux/osx" when running this function on a Windows operating system. 
#'@param parallelism.seed Either a numeric integer or NULL. This argument can be used to seed a random number generator which ensures reproducible calculations. This argument is effective regardless of whether parallel computations are being used.
#'@param n.optim A numeric integer. The number of times the Greedy Exchange Algorithm (Falk et al., 2014) is iterated to find an optimal design. Any integer greater than or equal to 1 is permissible, though larger values will produce more reliable results (while multiplying run-time). Defaults to 5.
#'@param n.draws A numeric integer. The number of Monte Carlo draws used to approximate the expected utility from the utility function per Muller (1999). Any values larger than 1 are permitted, though a minimum of 100 are recommended for the most stable results. 
#'@param extra.arguments A list of miscellaneous arguments and values which may be used to control the behaviour of the utility.function. See Details for more information.
#'@param ... Any additional arguments for the \code{foreach} iterator. The version of \code{foreach} from the package \code{doRNG} is used.
#'@return A list of four elements: 1) ssn.old, the original and unaltered ssn; 2) ssn.new, the original ssn modified such that it contains only the sites in the optimal design; 3) final.points, a vector of the locIDs or pids for the sampling sites present in the optimal design; and 4) utilities, a list of n.optim elements containing the expected utilities computed at every iteration of the Greedy Exchange Algorithm.
#' 
#'@details 
#' 
#' There are several aspects of this function that require explanations beyond the brief descriptions of the arguments above. 
#' \itemize{
#' \item n.points: this argument can be a single number or a vector. If the user supplies a single number, then, regardless of the number of isolated networks present in the ssn object, a total of n.points sites will be selected across ALL these networks. That is, the networks are not treated as separate design problems. However, if the user supplies a vector, things become more complicated. Firstly, the vector must be named. The names of the elements correspond to networks in the ssn argument and the element itself is the number of sites which should be selected within that network. The n.points argument may therefore look like this: \code{c("1" = 5, "3" = 6)}. In this example, the user is asking for 5 sites to be chosen from network 1 and 6 from network 3. It also shows that there is no need to select sites in every network. One final caveat is that no sites will appear in ssn.new for any networks which are skipped over.
#' \item utility.function: there are a number of pre-defined utility functions such as \code{DOptimality} and \code{KOptimality}. See THIS PAGE for an exhaustive list. Instructions for creating user-defined utility functions can also be found through that link.
#' \item prior.parameters: the number of elements in the prior.parameters list must be the same as the number of covariance parameters in the fitted \code{glmssn} which is also passed to this function. For example, if the only covariance parameter is the nugget (i.e. there is no spatial autocorrelation), then prior.parameters would be a list with one element. Each list element must be structured as an anonymous function as follows: \code{prior.parameters[[i]] <- function(x) runif(x)}. The random sample function does not have to be \code{runif}. The key message is that the function must only have a single argument, which is the number of draws to be taken from a random sample function. The function \code{\link{constructLogNormalCovPriors}} is able to construct lists of log-normal priors from \code{glmssn} objects. 
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
#' p <- constructLogNormalCovPriors(m1)
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
  utility.function,
  prior.parameters,
  n.cores = 1,
  parallelism = "none",
  parallelism.seed = NULL,
  n.optim = 5, 
  n.draws = 500,
  extra.arguments = list(),
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
  if(n.draws < 2){
    stop("The argument n.draws must be a positive integer greater than 1. Please make sure this is the case. Note that the recommended minimum value for n.draws is 100.")
  }
  if(!(is.null(parallelism.seed) | is.numeric(parallelism.seed))){
    stop("The argument parallelism.seed should either be left as NULL or changed to a positive integer.")
  }
  if(parallelism.seed <= 0){
    stop("The argument parallelism.seed must be a POSITIVE integer.")
  }
  if(floor(parallelism.seed) != parallelism.seed){
    stop("The argument parallelism.seed must be an INTEGER such that floor(parallelism.seed) is equal to parallelism.seed")
  }
  
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
  by.locID <- length(unique(ssn.obs.data$locID)) != length(ssn.obs.data$pid) 
  
  # If locIDs, map locIDs to PIDs
  if(by.locID){
    message("Multiple pids per locID have been detected. The resulting design will be for points with unique locIDs.")
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
  obs.X  <- model.matrix(m.form, ssn.obs.data)
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
    for(j in 1:n.optim){
      
      # Take random sample of points
      r.point <- sample(c.point, n.final, FALSE)
      if(by.locID){
        r.point.eval <- c()
        for(z in 1:n.final){
          ind <- match(r.point[z], lIDs)
          r.point.eval <- c(r.point.eval, locIDs[[ind]])
        }
      } else {
        r.point.eval <- r.point
      }
      
      # Evaluate utility once, store information
      u <- u.star <- utility.function(ssn, glmssn, r.point.eval, prior.parameters, n.draws, extra.arguments)
      r <- r.point
      
      # Initialise loop
      designs <- list(r.point)
      d.star <- designs[[1]]
      u.all <- c(u)
      condition <- TRUE
      while(condition){
        
        # Construct list of designs to evaluate
        u.this.run <- c()
        designs.this.run <- designs.this.run.eval <- list()
        for(k in 1:n.final){
          
          ind <- !(c.point %in% r.point) 
          not.in <- c.point[ind]
          
          # Replace k'th position with by cycling through remaining values
          for(l in 1:length(not.in)){
            
            r.point[k] <- not.in[l]
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
            
            designs.this.run <- append(designs.this.run, list(r.point))
            designs.this.run.eval <- append(designs.this.run.eval, list(r.point.eval))
            
          }
          
        }
        
        # Evaluate utilities HERE in parallel
        n.eval <- length(designs.this.run.eval)
        u.this.run <- foreach(i = 1:n.eval, .combine = c, .packages = c("SSN", "SSNdesign"), ...) %dopar% {
          utility.function(ssn, glmssn, designs.this.run.eval[[i]], prior.parameters, n.draws, extra.arguments)
        }
        
        u.all <- c(u.all, u.this.run)
        designs <- append(designs, designs.this.run)
        
        # Check logical condition
        condition <- any(u.this.run > u.star)
        if(condition){
          u.star <- max(u.this.run)[1]
          ind <- which(u.this.run == u.star)
          d.star <- designs[[ind[1]]]
          r.point <- d.star
        } 
        
      }
      
      # Store best
      best.u.in.opt <- c(best.u.in.opt, u.star)
      best.in.optim <- append(best.in.optim, list(c(d.star, x.point)))
      us.per.optims <- append(us.per.optims, list(u.all))
      
    }
    
    ind <- which(best.u.in.opt == max(best.u.in.opt))
    f.p.nets[[i]] <- best.in.optim[[ind[1]]]
    f.u.nets[[i]] <- us.per.optims
    
  }
  
  # Collapse results into one list
  f.p <<- unlist(f.p.nets)
  
  # Subset based on list
  if(by.locID){
    ssn.new <- subsetSSN(ssn, new.ssn.path, locID %in% f.p)
  } else {
    ssn.new <- subsetSSN(ssn, new.ssn.path, pid %in% f.p)
  }
  
  # Return list, subset, and diagnostics
  
  return(list(ssn.old = ssn, ssn.new = ssn.new, final.points = f.p, utilities = f.u.nets))
  
}