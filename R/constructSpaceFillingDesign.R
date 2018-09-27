constructSpaceFillingDesign <- function(ssn, new.ssn.path, n.points, n.optim = 3, type = "maximin", euclidean.distance = FALSE, p = 10){
  
  # Early exits 
  if(any(n.points == 1)){
    stop("It does not make sense to construct a space-filling design with only one point.")
  }
  if(!(type %in% c("maximin", "morris.mitchell"))){
    stop("The argument type must be one of maximin or morris.mitchell")
  }
  
  # Detect whether we're dealing with PIDs or locIDs
  by.locID <- length(unique(getSSNdata.frame(ssn)$locID)) != length(getSSNdata.frame(ssn)$pid) 
  
  # If locIDs, map locIDs to PIDs
  if(by.locID){
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
            
            # Evaluate utility
            if(type == "maximin"){
              u.kl <- spaceFillingMaxiMin(ssn, NULL, r.point.eval, NULL, 0, extra.arguments)
            } else{
              u.kl <- spaceFillingMorrisMitchell(ssn, NULL, r.point.eval, NULL, 0, extra.arguments)
            }
            u.this.run <- c(u.this.run, u.kl)
            
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
    f.p.nets[[i]] <- best.in.optim[[ind]]
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