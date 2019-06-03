#' @inherit DOptimality
#'@export
EKOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  # Copy ssn
  ssn2 <- ssn

  # Cut down SSN to contain only the design and prediction points
  ind.x <- ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% design.points
  ind.c <- row.names(ssn@obspoints@SSNPoints[[1]]@point.coords) %in% design.points
  ind.n <- row.names(ssn@obspoints@SSNPoints[[1]]@network.point.coords) %in% design.points
  ssn2@obspoints@SSNPoints[[1]]@network.point.coords <- ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind.n, ]
  ssn2@obspoints@SSNPoints[[1]]@point.coords <- ssn@obspoints@SSNPoints[[1]]@point.coords[ind.c, ]
  ssn2@obspoints@SSNPoints[[1]]@point.data <- ssn@obspoints@SSNPoints[[1]]@point.data[ind.x, ] 
  
  # Cut down SSN pred points similarly
  #indp.x <- ssn@predpoints@SSNPoints[[1]]@point.data$pid %in% row.names(extra.arguments$Matrices.prd$d)
  indp.c <- row.names(ssn@predpoints@SSNPoints[[1]]@point.coords) %in% row.names(extra.arguments$Matrices.prd$d)
  cds.prd <- ssn2@predpoints@SSNPoints[[1]]@point.coords[indp.c, ]
  colnames(cds.prd) <- c("x", "y")
  
  # Get other model parameters
  td <- glmssn$args$useTailDownWeight
  cm <- glmssn$args$CorModels
  un <- glmssn$args$use.nugget
  ua <- glmssn$args$use.anisotropy
  re <- glmssn$sampInfo$REs
  
  # Cut down matrices involving the observations
  mat <- extra.arguments$Matrices.Obs
  ind.mat <- row.names(mat$d) %in% design.points
  mat$d <-  mat$d[ind.mat, ind.mat]
  mat$a <-  mat$a[ind.mat, ind.mat]
  mat$b <-  mat$b[ind.mat, ind.mat] 
  mat$w <-  mat$w[ind.mat, ind.mat]
  extra.arguments$Matrices.obs <- mat
  net.zero.obs <- extra.arguments$net.zero.obs[ind.mat, ind.mat]
  
  # Do the same for the pxo matrix
  mat <- extra.arguments$Matrices.pxo
  mat$d <-  mat$d[ind.mat, ]
  mat$a <-  mat$a[ind.mat, ]
  mat$b <-  mat$b[ind.mat, ] 
  mat$w <-  mat$w[ind.mat, ]
  extra.arguments$Matrices.pxo <- mat
  net.zero.pxo <- extra.arguments$net.zero.pxo[ind.mat, ]
  
  # Get simulated empirical fixed effects to generate prior predictive datasets
  fep <- extra.arguments$Empirical.FEP
  
  # Estimate utility
  EK <- vector("numeric", n.draws)
  mod.formula <- as.character(glmssn$args$formula)
  nterms <- length(mod.formula)
  mod.formula <- as.formula(mod.formula[c(1, 3:nterms)])
  for(i in 1:n.draws){
    # Simulate data from simulated FE and CovParms values
    ssn.i <- SimulateOnSSN_minimal(
      ssn.object = ssn2,
      ObsSimDF = ssn2@obspoints@SSNPoints[[1]]@point.data,
      PredSimDF = NULL,
      PredID = NULL,
      formula = mod.formula,
      coefficients = fep[i, ],
      CorModels = glmssn$args$CorModels,
      use.nugget = glmssn$args$use.nugget,
      CorParms = prior.parameters[i, ],
      use.anisotropy = glmssn$args$use.anisotropy,
      addfunccol = glmssn$args$addfunccol,
      useTailDownWeight = glmssn$args$useTailDownWeight,
      family = glmssn$args$family,
      matrices.obs = extra.arguments$Matrices.obs,
      net.zero.obs = net.zero.obs,
      matrices.preds = extra.arguments$Matrices.prd,
      net.zero.preds = extra.arguments$net.zero.prd,
      matrices.predsxobs = extra.arguments$Matrices.pxo,
      net.zero.predsxobs = net.zero.pxo
    )$ssn.object
    # Fit model to simulated data
    mdl.tmp <- glmssn_minimal(
      formula = glmssn$args$formula,
      ssn.object = ssn.i,
      family = glmssn$args$family,
      CorModels = glmssn$args$CorModels,
      use.nugget = glmssn$args$use.nugget,
      use.anisotropy = glmssn$args$use.anisotropy,
      addfunccol = glmssn$args$addfunccol,
      useTailDownWeight = glmssn$args$useTailDownWeight,
      d = extra.arguments$Matrices.obs$d,
      a = extra.arguments$Matrices.obs$a,
      b = extra.arguments$Matrices.obs$b,
      c = extra.arguments$Matrices.obs$c,
      w = extra.arguments$Matrices.obs$w,
      n = net.zero.obs
    )
    # Obtain kriging variances at prediction sites from this model
    n.preds.uncertainty <- ncol(mdl.tmp$ssn.object@predpoints@SSNPoints[[1]]@point.data) + 1 # based on behaviour of predict.glmssn
    EK.i <- predict.glmssn_minimal(
      mdl.tmp,
      extra.arguments$Matrices.obs$a,
      extra.arguments$Matrices.obs$b,
      extra.arguments$Matrices.obs$w,
      net.zero.obs,
      extra.arguments$Matrices.obs$d,
      ssn2@obspoints@SSNPoints[[1]]@point.coords,
      extra.arguments$Matrices.prd$a,
      extra.arguments$Matrices.prd$b,
      extra.arguments$Matrices.prd$w,
      net.zero.prd,
      extra.arguments$Matrices.prd$d,
      cds.prd,
      extra.arguments$Matrices.pxo$a,
      extra.arguments$Matrices.pxo$b,
      extra.arguments$Matrices.pxo$w,
      net.zero.pxo,
      extra.arguments$Matrices.pxo$d
    )$ssn.object@predpoints@SSNPoints[[1]]@point.data[,n.preds.uncertainty]^2
    # Sum, invert
    EK[i] <- 1/sum(EK.i[!is.nan(EK.i)])
  }
  if(any(is.nan(EK))){
    return(c(-1e9))
  } else {
    return(mean(EK))
  }
}
