subsetByPID <- function(ssn, path, obs, preds){
  
  ssn <- subsetSSN(ssn, path, subset = ssn@obspoints@SSNPoints[[1]]@point.data$pid %in% obs)
  if(!is.null(preds)){
    ssn <- subsetPreds(ssn, pid %in% preds)
    createDistMat(ssn, "preds", T, T)
  } else {
    createDistMat(ssn, o.write = TRUE)
  }
  return(ssn)
  
}