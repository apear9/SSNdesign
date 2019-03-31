#' A function to construct a list of independent priors on covariance parameters from a glmssn object
#' 
#' @description 
#' 
#' The main design function in \code{SSNdesign} is called \code{optimiseSSNDesign}, and it has an argument specifying a list of independent priors on the covariance parameters of a fitted glmssn object.
#' 
#' @param glmssn A fitted glmssn object, preferably with \code{optimOutput$hessian}.
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
#' s <- SimulateOnSSN(s, getSSNdata.frame(s), formula = ~ 1, CorModels = c("Spherical.tailup"), CorParms = c(1, 2, 0.1), addfunccol = "addfunccol")$ssn.object
#' 
#' # Fit model
#' m <- glmssn(s, formula = Sim_Values ~ 1, CorModels = c("Spherical.tailup"), addfunccol = "addfunccol")
#' 
#' # Construct log-normal priors
#' p <- constructLogNormalPriors(m)
#' 
#' @export
constructLogNormalPriors <- function(glmssn, std){
  
  # Check inputs
  if(class(glmssn) != "glmssn"){
    stop("The argument glmssn must be an object of class glmssn.")
  }
  if(missing(std)){
    # Can either use observed or expected fisher info
    if(is.null(glmssn$optimOutput$hessian)){
      warning("The argument glmssn does not contain glmssn$optimOutput$hessian, so the expected fisher information will be used instead.")
      # theta <- glmssn$estimates$theta[1:length(glmssn$estimates$theta)]
      std <- theta.se <- getExpectedSE(glmssn, log.scale = TRUE)
      # std <- log.theta.se <- (1/theta) * theta.se
    } else {
      # Use delta method to find the se of the log(thetas)
      # inv.hessian <- solve(glmssn$optimOutput$hessian)
      # theta <- glmssn$estimates$theta[1:length(glmssn$estimates$theta)]
      # theta.se <- sqrt(diag(inv.hessian))
      std <- theta.se <- getCalculatedSE(glmssn, log.scale = TRUE)
    }
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