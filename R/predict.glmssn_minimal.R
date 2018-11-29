predict.glmssn <- function(object, a.mat.obs, b.mat.obs, w.mat.obs, n.zero.obs, d.hydro.obs, cds.obs,
                           a.mat.prd, b.mat.prd, w.mat.prd, n.zero.prd, d.hydro.prd, cds.prd, yeet){
  if(length(object$estimates$Warnlog) > 0 &&
     length(grep("Algorithm diverging",object$estimates$Warnlog)) > 0)
    stop("No predictions for diverging algorithm")
  datao <- object$ssn.object@obspoints@SSNPoints[[1]]@point.data
  ocoord <- object$ssn.object@obspoints@SSNPoints[[1]]@point.coords
  theta <- object$estimates$theta
  ind <- object$sampinfo$ind.obs
  nobs <- length(ind)
  CorModels <- object$args$CorModels
  useTailDownWeight <-object$args$useTailDownWeight
  REs <- object$sampinfo$REs
  distord <- order(as.integer(as.character(datao[,"netID"])),
                   datao[,"pid"])
  
  a.mat <- a.mat
  b.mat <- b.mat
  net.zero <- n.zero
  w.matrix <- w.mat
  dist.hydro <- d.hydro
  xcoord <- NULL
  ycoord <- NULL
  xyobs <- NULL
  xypred <- NULL
  rnames <- NULL
  
  # create Euclidean matrix among observed data
  xyobs <- ocoord
  x.samp <- ocoord[distord,1,drop = F]
  y.samp <- ocoord[distord,2,drop = F]
  xypred <- pcoord
  x.pred <- pcoord[, 1, drop = F]
  y.pred <- pcoord[, 2, drop = F]
  if(any(rownames(x.samp)!=rownames(a.mat)))
    stop("rownames of x.samp do not match rownames of a.mat")
  if(any(rownames(x.pred)!=colnames(a.mat)))
    stop("rownames of x.pred do not match colnames of a.mat")
  REPs <- NULL
  if(!is.null(REs)) {
    REnames <- names(REs)
    for(ii in 1:length(REnames)) if(any(is.na(datap[,REnames[ii]])))
      stop("Cannot having missing values when creating random effects")
    REOs <- list()
    REPs <- list()
    ObsSimDF <- datao
    PredSimDF <- datap
    ## model matrix for a RE factor
    for(ii in 1:length(REnames)){
      #we'll add "o" to observed levels and "p" to prediction
      # levels so create all possible levels
      plevels <- unique(c(levels(PredSimDF[,REnames[[ii]]]),
                          paste("o",levels(ObsSimDF[,REnames[[ii]]]),sep = ""),
                          paste("p",levels(PredSimDF[,REnames[[ii]]]),sep = "")))
      # sites with prediction levels same as observation levels
      pino <- PredSimDF[,REnames[[ii]]] %in% ObsSimDF[,REnames[[ii]]]
      #add "o" to observed levels
      ObsSimDF[,REnames[[ii]]] <- paste("o",
                                        ObsSimDF[,REnames[[ii]]], sep = "")
      ObsSimDF[,REnames[[ii]]] <- as.factor(as.character(
        ObsSimDF[,REnames[[ii]]]))
      #add all possible levels to prediction data frame
      levels(PredSimDF[,REnames[[ii]]]) <- plevels
      # add "o" to prediction sites with observation levels
      if(any(pino)) PredSimDF[pino,REnames[[ii]]] <- paste("o",
                                                           PredSimDF[pino,REnames[[ii]]], sep = "")
      # add "p" to all predicition sites without observation levels
      if(any(!pino)) PredSimDF[!pino,REnames[[ii]]] <- paste("p",
                                                             PredSimDF[!pino,REnames[[ii]]], sep = "")
      PredSimDF[,REnames[[ii]]] <- as.factor(as.character(
        PredSimDF[,REnames[[ii]]]))
      # now get down to just levels with "o" & "p" added
      blevels <- unique(c(levels(ObsSimDF[,REnames[[ii]]]),
                          levels(PredSimDF[,REnames[[ii]]])))
      ObsSimDF[,REnames[[ii]]] <- factor(ObsSimDF[,REnames[[ii]]],
                                         levels = blevels, ordered = FALSE)
      PredSimDF[,REnames[[ii]]] <- factor(PredSimDF[,REnames[[ii]]],
                                          levels = blevels, ordered = FALSE)
      # now ordering of factors in Z matrices should be compatible
      # with obs x obs Z matrices
      REOs[[ii]] <- model.matrix(~ObsSimDF[distord,
                                           REnames[[ii]]] - 1)
      REPs[[ii]] <- model.matrix(~PredSimDF[,
                                            REnames[[ii]]] - 1)
      rownames(REOs[[ii]]) <- datao[distord,"pid"]
      rownames(REPs[[ii]]) <- datap[,"pid"]
      if(any(rownames(REOs[[ii]])!=rownames(a.mat)))
        stop("rownames RE for obs do not match rownames of a.mat")
      if(any(rownames(REPs[[ii]])!=colnames(a.mat)))
        stop("rownames RE for preds do not match colnames of a.mat")
    }
    ## corresponding block matrix
    for(ii in 1:length(REnames)) REPs[[ii]] <-
      REOs[[ii]] %*% t(REPs[[ii]])
  }
  Vpred <- makeCovMat(theta = theta, dist.hydro = dist.hydro,
                      a.mat = a.mat, b.mat = b.mat, w.matrix = w.matrix,
                      net.zero = net.zero, x.row = x.samp, y.row = y.samp,
                      x.col = x.pred, y.col = y.pred,
                      CorModels = CorModels,
                      useTailDownWeight = useTailDownWeight,
                      use.nugget = FALSE,
                      use.anisotropy = object$args$use.anisotropy, REs = REPs)
  Vpred <- Vpred[ind, , drop = F]
  
  # get a list of response and covariate names
  response.col <- object$args$zcol
  mod.names <- as.character(attr(terms(object$args$formula,
                                       data = object$ssn.object@obspoints@SSNPoints[[1]]@point.data),"variables"))
  # get the number of names ( + 1, as the first is always "list")
  nc.tmp <- length(mod.names)
  # if there are any covariates ...
  ind.allcov <- rep(TRUE, times = length(datap[,1]))
  if(nc.tmp > 2) {
    # create a FALSE for a record with missing values of the covariates
    for(i in 3:nc.tmp) ind.allcov <- ind.allcov & !is.na(datap[,mod.names[i]])
  }
  # prediction sample size without missing covariates
  np.allcov <- sum(ind.allcov)
  # add response variable as -1 in datap
  datap[,response.col] <- NA
  # add prediction standard errors as NAs in datap
  datap[,paste(response.col,".predSE", sep = "")] <- NA
  # remove records that had any missing values for any of the covariates
  datap1 <- datap[ind.allcov,]
  # add response variable as -1 in datap
  datap1[,response.col] <- -1
  # add prediction standard errors as NAs in datap
  datap1[,paste(response.col,".predSE", sep = "")] <- NA
  
  formula <- object$args$formula
  # create design matrix for prediction data set
  mf <- model.frame(formula, data = datap1)
  mt <- attr(mf, "terms")
  Xpred <- model.matrix(mt, mf, contrasts)
  Xpred <- Xpred[,object$sampinfo$cutX1toX2, drop = F]
  
  # get the sum of partial sills
  sumparsil <- sum(theta[attr(theta,"type") == "parsill"])
  Vi <- object$estimates$Vi
  covb <- object$estimates$covb
  Xobs <- object$sampinfo$X
  z <- object$sampinfo$z
  n <- object$sampinfo$obs.sample.size
  p <- object$sampinfo$rankX
  
  parsilvec <- rep(sumparsil, times = length(Vpred[1,]))
  
  if(object$args$family == "poisson") {
    beta.hat <- object$estimates$betahat
    eta.hatp <- Xpred[ind.allcov,] %*% beta.hat
    eta.hato <- Xobs %*% beta.hat
    #diagonal elements of Delta~^{-1} of my manuscript
    Del.ip <- as.vector(1/exp(eta.hatp))
    Del.io <- as.vector(1/exp(eta.hato))
    #diagonal elements of A^(1/2) of my manuscript
    A.5p <- as.vector(sqrt(exp(eta.hatp)))
    A.5o <- as.vector(sqrt(exp(eta.hato)))
    Vpred <- t((Del.ip*A.5p)*t(Del.io*A.5o*Vpred))
    parsilvec <- sumparsil*(A.5p*Del.ip)^2
  }
  if(object$args$family == "binomial") {
    beta.hat <- object$estimates$betahat
    eta.hatp <- Xpred[ind.allcov,] %*% beta.hat
    eta.hato <- Xobs %*% beta.hat
    #diagonal elements of Delta~^{-1} of my manuscript
    Del.ip <- as.vector((1 + exp(eta.hatp))^2/exp(eta.hatp))
    Del.io <- as.vector((1 + exp(eta.hato))^2/exp(eta.hato))
    #diagonal elements of A^(1/2) of my manuscript
    A.5p <- as.vector(sqrt(exp(eta.hatp)/(1 +
                                            exp(eta.hatp))^2))
    A.5o <- as.vector(sqrt(exp(eta.hato)/(1 +
                                            exp(eta.hato))^2/object$sampinfo$trialsvec))
    Vpred <- t((Del.ip*A.5p)*t(Del.io*A.5o*Vpred))
    parsilvec <- sumparsil*(A.5p*Del.ip)^2
  }
  
  M <- rbind(Vpred, t(Xpred), parsilvec)
  XXSiXi <- Xobs %*% covb
  XSi <- t(Xobs) %*% Vi
  pred.out <- t(apply(M, 2, UK4Apply, covb = covb,
                      XXSiXi = XXSiXi, XSi = XSi, Vi = Vi, z = z, n = n, p = p))
  datap1[,response.col] <- pred.out[,1]
  datap1[,paste(response.col,".predSE", sep = "")] <- pred.out[,2]
  datap[ind.allcov,] <- datap1
  #put the predictions in the predicted data data.frame and return the SSN object
  for(i in 1:length((object$ssn.object@predpoints@SSNPoints)))
    if(object$ssn.object@predpoints@ID[i] == predpointsID){
      object$ssn.object@predpoints@SSNPoints[[i]]@point.data <-
        datap[order(distordp),]
    }
  #object$args$predpointsID <- predpointsID
  #class(object) <- "glmssn.predict"
  
  object
  
}

