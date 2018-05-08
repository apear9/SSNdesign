#' Pull the design diagnostics into a dataframe.
#'
#'@description 
#'
#' Creates a dataframe of utility values across all networks for all exchange algorithm hyper-iterations from the list output by the main design functions.  
#' 
#'@usage 
#' 
#' \code{ssndesign2DF(ssndesign, summarise = TRUE)}
#' 
#'@param ssndesign The list object output by one of the main design functions.
#'@return A data.frame object.
#'
#'@details 
#'
#'This function does stuff.
#'  
#'@export
ssndesign2DF <- function(ssndesign, summarise = TRUE){
  
  # Check the input
  if(!is.list(ssndesign)){
    stop("The argument ssndesign must be the list output by one of the main design functions.")
  }
  
  # Find the number of networks to plot traces for
  n.net <- length(ssndesign$utility.values) # List has one element per network
  
  # Find number of hyper iterations
  K <- length(ssndesign$utility.values[[1]]) # Each sublist contains K lists; one for each kth hyper-iteration.
  
  # Find whether networks treated separately
  
  sept <- ssndesign$networks.separate
  
  # Then find the number of iterations for the for loops in each iterations
  
  n.pts.orig <- ndpoints(ssndesign$ssn.object.old)
  n.pts.new  <- ndpoints(ssndesign$ssn.object.new)
  
  if(sept){
    num.evals.per.iter <- n.pts.new$byNetwork * (n.pts.orig$byNetwork - n.pts.new$byNetwork)
  } else {
    num.evals.per.iter <- n.pts.new$inTotal * (n.pts.orig$inTotal - n.pts.new$inTotal)
  }
  
  # Pull into dataframe
  df <- data.frame()
  
  # Two cases: summarise = TRUE -> summary stats. summarise = FALSE -> raw data
  if(summarise){
    
    for(i in 1:n.net){
      
      for(j in 1:K){
        
        # Few things to find here
        # 1. The random start for this kth hyper-iteration
        # 2. The number of times the while loop ran
        n.j <- length(ssndesign$utility.values[[i]][[j]]) - 1 # The first value is always a random design
        n.iter <- n.j/num.evals.per.iter[i]
        
        # Begin constructing dataframe
        df.tmp <- data.frame(
          min = rep(0, n.iter + 1),
          median = rep(0, n.iter + 1),
          max = rep(0, n.iter + 1),
          K = rep(j, n.iter + 1),
          itr = 0:n.iter, # Zero is the random start
          network = rep(i, n.iter + 1)
        )
        
        # Find mins, maxes per iteration
        df.tmp[1, 1:3] <- ssndesign$utility.values[[i]][[j]][1]
        
        for(k in 2:(n.iter + 1)){
          indices.k <- (k + (k - 2) * num.evals.per.iter[i]):((k - 1) * num.evals.per.iter[i])
          df.tmp[k, 1] <- min(ssndesign$utility.values[[i]][[j]][indices.k])
          df.tmp[k, 2] <- median(ssndesign$utility.values[[i]][[j]][indices.k])
          df.tmp[k, 3] <- max(ssndesign$utility.values[[i]][[j]][indices.k])
        }
        
        df <- rbind(df, df.tmp)
        
      }
      
    }
    
  } else {
    
    for(i in 1:n.net){
      
      for(j in 1:K){
        
        # Few things to find here
        # 1. The random start for this kth hyper-iteration
        # 2. The number of times the while loop ran
        n.j <- length(ssndesign$utility.values[[i]][[j]]) - 1 # The first value is always a random design
        n.iter <- n.j/num.evals.per.iter[i]
        
        # Begin constructing dataframe
        df.tmp <- data.frame(
          utility = ssndesign$utility.values[[i]][[j]],
          K = rep(j, n.j + 1),
          itr = c(0, rep(1:n.iter, each = num.evals.per.iter[i])), 
          network = rep(i, n.j + 1)
          )
          
          df <- rbind(df, df.tmp)
          
      }
      
    }
    
  }
  
  ## Spit out the resulting dataframe
  
  return(df)
  
}

