design.validation.predictions <- function(full.ssn, designs.list, glmssn, dist.type, n.sims){
  
  ## Check that prediction sites are present
  any.preds <- unlist(
    lapply(
      designs.list,
      anyPreds
    )
  )
  if(!all(any.preds)){
    stop("No prediction points are present in the SpatialStreamNetwork.")
  }
  
  ## Find number of designs to validate
  
  n.designs <- length(designs.list)
  
  ## Create formula for model-fitting
  form.old <- as.character(glmssn$args$formula)
  form.new <- paste("SIM", form.old[3], sep = "~")
  form.new <- as.formula(form.new)
  
  ## Do simulations
  results <- matrix(0, nrow = n.designs, ncol = n.sims)
  for(i in 1:n.sims){
    
    # Simulations
    ssn.sim.i <- simulateFromSSNM(full.ssn, glmssn)
    sim.obs <- getSSNdata.frame(ssn.sim.i)
    sim.prd <- getSSNdata.frame(ssn.sim.i, "preds")
    
    ## Do validation for this set of simulated data
    results.sim.i <- vector("numeric", n.designs)
    for(j in 1:n.designs){
      
      # Extract simulated values
      rn.obs <- row.names(sim.obs)
      rn.prd <- row.names(sim.prd)
      pid.obs <- designs.list[[j]]@obspoints@SSNPoints[[1]]@point.data$pid
      pid.prd <- designs.list[[j]]@predpoints@SSNPoints[[1]]@point.data$pid
      ind.obs <- match(pid.obs, rn.obs)
      ind.prd <- match(pid.prd, rn.prd)
      designs.list[[j]]@obspoints@SSNPoints[[1]]@point.data$SIM <- sim.obs$Sim_Values[ind.obs]
      designs.list[[j]]@predpoints@SSNPoints[[1]]@point.data$SIM <- sim.prd$Sim_Values[ind.prd]
      
      # Fit model
      results.sim.i[j] <- tryCatch(
        expr = {
          model <- glmssn(
            formula = form.new,
            ssn.object = designs.list[[j]],
            family = glmssn$args$family,
            use.nugget = glmssn$args$use.nugget,
            CorModels = glmssn$args$CorModels,
            use.anisotropy = glmssn$args$use.anisotropy,
            addfunccol = glmssn$args$addfunccol,
            trialscol = glmssn$args$trialscol,
            EstMeth = glmssn$args$EstMeth,
            useTailDownWeight = glmssn$args$useTailDownWeight,
            trans.power = glmssn$args$trans.power,
            trans.shift = glmssn$args$trans.shift,
            control = glmssn$args$control
          )
          # Predict from model
          preds <- predict.glmssn(model, "preds")$ssn.object
          
          # Calculate distance between simulated and predicted values
          preds.true <- getSSNdata.frame(designs.list[[j]], "preds")
          preds.pred <- getSSNdata.frame(preds, "preds")
          sim <- preds.true$SIM
          pred <- preds.pred$SIM
          results.sim.i[j] <- dist.type(sim, pred)
          
        }, 
        
        error = function(e){NA} 
      )
      
    }
    
    ## Store results
    results[,i] <- results.sim.i
    
  }
  
  ## Construct and return list of results
  
  list(
    Raw = results,
    Failures = rowSums(is.na(results)),
    Avg = rowMeans(results, na.rm = TRUE)
  )
  
}
