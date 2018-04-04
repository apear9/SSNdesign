#' A function for fitting SpatialStreamNetwork models inside utility functions
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

