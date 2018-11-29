#' A function to construct a list of independent priors on covariance parameters from a glmssn object
#' 
#' @description 
#' 
#' The main design function in SSNdesign, \code{optimiseSSNDesign}, has an argument specifying a list of independent priors on the covariance parameters of a fitted glmssn object.
#' 
#' @usage 
#' 
#' \code{constructLogNormalCovPriors(glmssn, std = c(0.25, 0.5))}
#'
#' @param glmssn A fitted glmssn object.
#' @param std A numeric vector. This vector specifies the standard deviation of the log-normal priors. See Details for more information.
#' @return A list whose elements are functions parameterised in terms of x, the number of Monte Carlo draws to be taken from the priors when evaluating the expected utility. 
#' 
#' @details 
#' 
#' The argument std can be of length 1, in which case all covariance parameters (parsill, range) will have the same value for standard deviation. Alternatively, it can be of length two, in which case the first element will be used for the partial sill parameters and the second for range parameters. Alternatively, it can be of the same length as glmssn$estimates$theta in which the elements will be matched to their corresponding parameters.
#' 
#' @examples 
#' 
#' # Create SSN
#' s <- createSSN(10, binomialDesign(10), path = paste(tempdir(), "s.ssn", sep = "/"), importToR = TRUE)
#' 
#' # Simulate data
#' s <- SimulateOnSSN(s, getSSNdata.frame(s), formula = ~ 1, CorModels = c("Exponential.tailup"), CorParms = c(1, 2, 0.1), addfunccol = "addfunccol")$ssn.object
#' 
#' # Fit model
#' m <- glmssn(s, formula = Sim_Values ~ 1, CorModels = c("Exponential.tailup"), addfunccol = "addfunccol")
#' 
#' # Construct log-normal priors
#' p <- constructLogNormalCovPriors(m)
#' 
#' @export
constructLogNormalCovPriors <- function(glmssn, std = c(0.25, 0.5)){
  
  # Check inputs
  if(class(glmssn) != "glmssn"){
    stop("The argument glmssn must be an object of class glmssn.")
  }
  if(!is.numeric(std)){
    stop("The argument std must be numeric.")
  }
  
  # extract out estimates of covariance parameters
  cov <- glmssn$args$CorModels
  est <- log(glmssn$estimates$theta)
  n.cov <- length(cov)
  n.est <- length(est)
  n.std <- length(std)
  a.nug <- (glmssn$args$use.nugget) * 1
  
  # Check compatibility of arguments and construct true vector of standard deviations
  if(n.est < n.std){
    stop("Too many values have been provided for the argument std.")
  }
  if(n.std == 1){
    std.use <- rep(std, n.est)
  }
  if(n.std == 2){
    std.use <- rep(std, n.cov)
    if(a.nug){
      std.use[n.est] <- std[1]
    }
  }
  if(n.std == n.est){
    std.use <- std
  }
  
  # Loop to construct stuff by parsing
  functions <- vector("character", n.est)
  for(i in 1:n.est){
    functions[i] <- paste0("function(x){exp(rnorm(x,", est[i], ",", std.use[i],"))}")
  }
  functions <- paste(functions, collapse = ",")
  
  # Add final parts of command
  to.run <- paste0("list(", functions, ")")
  
  # Construct list by parsing command
  priors <- eval(
    parse(
      text = to.run
    )
  )
  
  return(priors)
  
}