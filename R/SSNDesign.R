#' SSNdesign: An R package for optimal and adaptive sampling designs on stream networks
#'
#' The SSNdesign package builds on the SSN package to provide a framework for solving design problems on SpatialStreamNetwork objects.
#'
#' @section SSNdesign functions:
#' the SSNdesign package contains the following functions:
#' - constructTotalMatrix
#' - CPOptimality
#' - doAdaptiveDesign
#' - doAdaptiveDesignParallel
#' - DOptimality
#' - EDOptimality
#' - EKOptimality
#' - findOptimalDesign
#' - findOptimalDesignParallel
#' - generateSites
#' - getBIDtables
#' - getImportantMatrices.obs
#' - getImportantMatrices.pxo
#' - getStreamDistMatsInOrder.predsxobs
#' - getStreamDistMatsInOrder
#' - glmssn_minimal
#' - hardCoreDesign
#' - importStreams
#' - KOptimality
#' - ndpoints
#' - nnetwork
#' - nppoints
#' - searchAllDesigns
#' - sequentialCPOptimality
#' - sequentialDOptimality
#' - sequentialEDOptimality
#' - SimulateOnSSN_minimal
#' - stableInverse
#' - systematicDesign
#'
#' @docType package
#' @name SSNdesign
#' @import foreach
#' @import iterators
#' @import itertools
#' @import doParallel
#' @import SSN
#' @import DBI
#' @import shp2graph
#' @import maptools
#' @import igraph
#' @import rgdal
#' @import RSQLite
#' @import stringr
#' @import doRNG
#' @importFrom dplyr left_join inner_join
#' @importFrom ggplot2 fortify
#' @import spsurvey
#' @importFrom MASS mvrnorm
#' @importFrom graphics plot lines
#' @importFrom methods is new
#' @importFrom stats as.formula contrasts dist glm model.frame model.matrix optim optimize rbinom rnorm rpois sd terms
#' @importFrom utils read.table
#' @importFrom sp spTransform SpatialPointsDataFrame SpatialLines SpatialLinesDataFrame
#' @importFrom parallel makeCluster stopCluster
NULL
## NULL
