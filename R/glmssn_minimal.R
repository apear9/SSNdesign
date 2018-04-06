#' A function for fitting SpatialStreamNetwork models inside utility functions
#' 
#'@description
#'
#'This function is the equivalent of glmssn from the package SSN. However, unlike that function, this function requires the user to input the distance and weights matrices for the observed sites over which the model is being fit.
#'
#'@usage
#'
#'\code{glmssn_minimal(..., d, a, b, c, w)}
#'
#'@param ... Arguments for the function \code{glmssn} from the package SSN.
#'@param d,a,b,c,w Distance matrices (d, a, b), a connectivity matrix (c), and a weights matrix (w).
#'@return An object of class glmssn, as described for the function \code{glmssn} in the package SSN.
#'
#'@details
#'
#' This function is primarily intended for use inside user-defined utility functions, when the utility function requires that models be fit to the data for the design points under consideration. The way this function is set up avoids doubling up on the computations required to generate the distance matrices (which will already have been created prior to the utility function being called) and errors associated with attempting to read in the distance matrices for singular networks.
#' 
#'@export 
glmssn_minimal <- function(formula, ssn.object,
                   family = "Gaussian",
                   CorModels = c("Exponential.tailup",
                                 "Exponential.taildown",
                                 "Exponential.Euclid"),
                   use.nugget = TRUE,
                   use.anisotropy = FALSE,
                   addfunccol = NULL,
                   trialscol = NULL,
                   EstMeth = "REML",	
                   useTailDownWeight = FALSE,
                   trans.power = NULL,
                   trans.shift = 0,
                   control = list(max.range.factor = 4, 
                                  trunc.pseudo = NULL, 
                                  maxiter.pseudo = 20,
                                  beta.converge = 1e-5),
                   d,
                   a,
                   b,
                   c,
                   w
)
{
  
  output <- glmssn1_minimal(formula=formula,
                    ssn.object=ssn.object,
                    family = family,
                    CorModels = CorModels,
                    useTailDownWeight = useTailDownWeight,
                    use.nugget = use.nugget,
                    use.anisotropy = use.anisotropy,
                    addfunccol = addfunccol,
                    trialscol = trialscol,
                    EstMeth = EstMeth,
                    trans.power = trans.power,
                    trans.shift = trans.shift,
                    control = control,
                    dist.hydro.data = d,
                    a.mat.data = a,
                    b.mat.data = b,
                    c.mat.data = c,
                    w.matrix.data = w
  )
  cl = match.call()
  output$args[["call"]] = cl
  return(output)
  
}

