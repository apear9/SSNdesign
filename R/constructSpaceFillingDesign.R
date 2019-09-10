#' Get a space-filling design for a Spatial Stream Network
#' 
#' @description 
#' 
#' This function finds a space-filling design on a Spatial Stream Network. It does so by optimising one of two space-filling utility functions: (1) the maximin utility function in \code{\link{spaceFillingMaxiMin}} or (2) the modified maximin utility proposed by Morris and Mitchell (1995) in \code{\link{spaceFillingMorrisMitchell}}. 
#' 
#' The Greedy Exchange Algorithm is used to find the optimal design.
#'  
#' @param ssn An object of class SpatialStreamNetwork. Its obspoints slot must contain all the potential sampling sites for the stream network.  
#' @param new.ssn.path A path to a folder where the result can be stored.
#' @param n.points A numeric or a named numeric vector specifying the size of the final design(s). See Details for more information.
#' @param type A character string. Either \code{"maximin"} or \code{"morris.mitchell"}. 
#' @param euclidean.distance A logical. This specifies whether Euclidean distance should be used to construct the space-filling designs, instead of hydrological distance within the stream network. Defaults to FALSE.
#' @param p A numeric. This is a weighting power used by the morris.micthell utility function. It is ignored when \code{type = "maximin"}. This number must be greater than or equal to 1. 
#' @param n.cores The number of CPUs which should be used when executing this function. This argument must agree with the argument parallelism. For example, if n.cores > 1 and parallelism is "none", this argument will be ignored and all computations will be performed sequentially. Defaults to 1. 
#' @param parallelism A character which must be one of "none", "windows", "linux/osx". These can be spelled with in any case or in any combination of cases. Note the argument must be selected appropriately for the operating system on the user's computer. Definitely do not select "linux/osx" when running this function on a Windows operating system. 
#' @param parallelism.seed Either a numeric integer or NULL. This argument can be used to seed a random number generator which ensures reproducible calculations. This argument is effective regardless of whether parallel computations are being used.
#' @param n.optim A numeric integer. The number of times the Greedy Exchange Algorithm (Falk et al., 2014) is iterated to find an optimal design. Any integer greater than or equal to 1 is permissible, though larger values will produce more reliable results (while multiplying run-time). Defaults to 5.
#' @param ... Any additional arguments for the \code{foreach} iterator. The version of \code{foreach} from the package \code{doRNG} is used.
#' @return A list of four elements: 1) ssn.old, the original and unaltered ssn; 2) ssn.new, the original ssn modified such that it contains only the sites in the optimal design; 3) final.points, a vector of the locIDs or pids for the sampling sites present in the optimal design; and 4) utilities, a list of n.optim elements containing the expected utilities computed at every iteration of the Greedy Exchange Algorithm.
#' 
#' @details 
#' 
#' The argument n.points can be a single number or a vector. If the user supplies a single number, then, regardless of the number of isolated networks present in the ssn object, a total of n.points sites will be selected across ALL these networks. That is, the networks are not treated as separate design problems. However, if the user supplies a vector, things become more complicated. Firstly, the vector must be named. The names of the elements correspond to networks in the ssn argument and the element itself is the number of sites which should be selected within that network. The n.points argument may therefore look like this: \code{c("1" = 5, "3" = 6)}. In this example, the user is asking for 5 sites to be chosen from network 1 and 6 from network 3. It also shows that there is no need to select sites in every network. One final caveat is that no sites will appear in ssn.new for any networks which are skipped over.
#' 
#' @references 
#' 
#' Falk, M., Pettitt, A., McGree, J.M. (2014). Sampling designs on stream networks using the pseudo-Bayesian approach. \emph{Environmental and Ecological Statistics}, \emph{21}(4), 751-773.
#' 
#' Morris, M.D. & Mitchell, T.J. (1995). Exploratory Designs for Computational Experiments. \emph{Journal of Statistical Planning and Inference}, \emph{43}, 381-402. 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#'  # Create a simulated SpatialStreamNetwork object
#'  s <- createSSN(50, binomialDesign(100), 
#'  path = tempPath("s.ssn"), importToR = TRUE)
#'  createDistMat(s)
#'  
#'  # Construct a space filling design
#'  space.filling <- constructSpaceFillingDesign(s, 
#'  tempPath("r.ssn"), 50, "maximin")
#'  
#'  # Plot result to check
#'  plot(space.filling$ssn.new)
#'  
#' }
#'   
#' @export
constructSpaceFillingDesign <- function(
  ssn, 
  new.ssn.path, 
  n.points, 
  type = "maximin", 
  euclidean.distance = FALSE, 
  p = 10,
  n.optim = 5, 
  n.cores = 1,
  parallelism = "none",
  parallelism.seed = NULL,
  ...
){
  
  # Early exits
  if(!isSSN(ssn)){
    stop("The argument SSN must be of class SpatialStreamNetwork.")
  }
  if(dir.exists(new.ssn.path)){
    stop("The SSN folder specified in the argument new.ssn.path already exists. Please choose a new one or delete the existing folder.")
  }
  if(any(n.points == 1)){
    stop("It does not make sense to construct a space-filling design with only one point.")
  }
  if(!(type %in% c("maximin", "morris.mitchell"))){
    stop("The argument type must be one of maximin or morris.mitchell")
  }
  if(!is.character(parallelism)){
    stop("The parallelism must be one of the following strings, with any capitalisation: none, windows, osx/linux.")
  }
  parallelism <- stringr::str_to_lower(parallelism)
  if(!(parallelism %in% c("none", "windows", "osx/linux"))){
    stop("The argument parallelism must be one of : none, windows, osx/linux.")
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
  
  # Detect whether we're dealing with PIDs or locIDs
  by.locID <- length(unique(getSSNdata.frame(ssn)$locID)) != length(getSSNdata.frame(ssn)$pid) 
  
  # If locIDs, map locIDs to PIDs
  if(by.locID){
    message("Multiple pids per locID have been detected. The resulting design will be for points with unique locIDs.")
    alls <- getSSNdata.frame(ssn)$locID
    lIDs <- unique(alls)
    pIDs <- getSSNdata.frame(ssn)$pid
    locIDs <- vector("list", length(lIDs))
    for(i in 1:length(lIDs)){
      locIDs[[i]] <- pIDs[alls == lIDs[i]]
    }
    names(locIDs) <- lIDs
  }
  
  # Find number of points in each network
  netIDs <- getSSNdata.frame(ssn)$netID
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
  
  # Extract out distance matrices and stuff
  if(euclidean.distance){
    d <- as.matrix(dist(ssn@obspoints@SSNPoints[[1]]@point.coords, diag = TRUE))
    rownames(d) <- colnames(d) <- getSSNdata.frame(ssn)$pid
  } else {
    d <- getStreamDistMatInOrder(ssn)
    d <- constructTotalMatrix(d)
    n.z <- d$net.zero
    d <- getImportantMatrices.obs(d$d.junc, rep(1, nrow(ssn@obspoints@SSNPoints[[1]]@point.data)))$d
    d <- d * n.z
  }
  extra.arguments <- list(Matrices.Obs = list(d = d), TND = upper.tri(d))
  if(type == "morris.mitchell"){
    extra.arguments$p <- p
  }
  
  # Construct sets of candidate points
  if(by.locID){
    if(length(n.points) != 1){
      c.points <- vector("list", length(n.points))
      for(i in 1:length(n.points)){
        ind <- getSSNdata.frame(ssn)$netID == d.nets.names[i]
        c.points[[i]] <- getSSNdata.frame(ssn)$locID[ind]
      }
    } else {
      c.points <- list(getSSNdata.frame(ssn)$locID)
    }
  } else {
    if(length(n.points) != 1){
      c.points <- vector("list", length(n.points))
      for(i in 1:length(n.points)){
        ind <- getSSNdata.frame(ssn)$netID == d.nets.names[i]
        c.points[[i]] <- getSSNdata.frame(ssn)$pid[ind]
      }
    } else {
      c.points <- list(getSSNdata.frame(ssn)$pid)
    }
  }
  
  # Loop through to optimise
  f.p.nets <- vector("list", length(n.points))
  f.u.nets <- vector("list", length(n.points))
  for(i in 1:length(n.points)){
    
    n.final <- n.points[i]
    c.point <- c.points[[i]]
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
      
      # Evaluate utility once
      if(type == "maximin"){
        u <- u.star <- spaceFillingMaxiMin(ssn, NULL, r.point.eval, NULL, 0, extra.arguments)
      } else {
        u <- u.star <- spaceFillingMorrisMitchell(ssn, NULL, r.point.eval, NULL, 0, extra.arguments)
      }
      d.star <- r.point
      
      # Initialise loop
      designs <- list(r.point)
      u.all <- c(u)
      condition <- TRUE
      while(condition){
        
        # Construct list of designs to evaluate
        u.this.run <- c()
        designs.this.run <- list()
        designs.this.run.eval <- list()
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
            } else {
              r.point.eval <- r.point
            }
            
            designs.this.run <- append(designs.this.run, list(r.point))
            designs.this.run.eval <- append(designs.this.run.eval, list(r.point.eval))
            
          }
          
        }
        
        # Evaluate utilities HERE in parallel
        n.eval <- length(designs.this.run.eval)
        if(type == "maximin"){
          u.this.run <- foreach(i = 1:n.eval, .combine = c, .packages = c("SSN", "SSNdesign"), ...) %dopar% {
            spaceFillingMaxiMin(ssn, NULL, designs.this.run.eval[[i]], NULL, 0, extra.arguments)
          }
        } else {
          u.this.run <- foreach(i = 1:n.eval, .combine = c, .packages = c("SSN", "SSNdesign"), ...) %dopar% {
            spaceFillingMorrisMitchell(ssn, NULL, designs.this.run.eval[[i]], NULL, 0, extra.arguments)
          }
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
      best.in.optim <- append(best.in.optim, list(d.star))
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