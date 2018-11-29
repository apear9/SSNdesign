design.validation.fixed.effects.estimation <- function(full.ssn, designs.list, glmssn, dist.type, n.sims){
  
  ## Get true values of standard errors in the model
  se_true <- diag(glmssn$estimates$covb)
  
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

    ## Do validation for this set of simulated data
    results.sim.i <- vector("numeric", n.designs)
    for(j in 1:n.designs){
      
      # Extract simulated values
      rn.obs <- row.names(sim.obs)
      pid.obs <- designs.list[[j]]@obspoints@SSNPoints[[1]]@point.data$pid
      ind.obs <- match(pid.obs, rn.obs)
      designs.list[[j]]@obspoints@SSNPoints[[1]]@point.data$SIM <- sim.obs$Sim_Values[ind.obs]
      
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
          
          # Extract se of FE out of model
          se_fit <- diag(model$estimates$covb) # Ask James whether diagonal is what I want
          
          # Calculate distance between fit and true models in terms of se of fe
          results.sim.i[j] <- dist.type(se_fit, se_true)
          
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