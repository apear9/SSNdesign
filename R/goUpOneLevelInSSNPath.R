goUpOneLevelInSSNPath <- function(path){
  
  # Get path info
  
  full <- file.path(path)
  abbr <- basename(full)
  
  # Up one level
  gsub(abbr, "", full)
  
}