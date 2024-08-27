#' @title UncertainSolver
#' @name UncertainSolver
#' @description solver to analyze sensitivity of parameters
#' @param ParentModule SBCore
#' @param vnamesDistSD  dataframe with columns vnames, distNames (see lhs package for possible distributions), secondPar
#' @param n samplesize 
#' @return States (i) (=mass)
#' @export
UncertainSolver = function(ParentModule, tol=1e-30) { 
  
  TheCore <- ParentModule$myCore
  
  
  
  
  
  TheCore$mutateVars(Updated)
  
  #update core and solve
  TheCore$UpdateDirty(uniqvNames)
  ParentModule$PrepKaasM()
  ParentModule$PrepemisV()
  
  
  
  
  
  SB.K = ParentModule$SB.k
  vEmis = ParentModule$emissions
  
  solve(SB.K, -vEmis, tol = tol)
  
}