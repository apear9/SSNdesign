#' @inherit DOptimality
#' @export
CPDOptimality <- function(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments){
  
  CPOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments) + DOptimality(ssn, glmssn, design.points, prior.parameters, n.draws, extra.arguments)
  
}