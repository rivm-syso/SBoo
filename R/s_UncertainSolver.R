#' @title UncertainSolver
#' @name UncertainSolver
#' @description solver to analyze sensitivity of parameters
#' @param ParentModule SBCore
#' @param vnamesDistSD  dataframe with columns vnames, distNames (see lhs package for possible distributions), secondPar
#' @param n samplesize 
#' @return States (i) (=mass)
#' @export
UncertainSolver = function(ParentModule, tol=1e-30) { 
  browser()
  TheCore <- ParentModule$myCore
  sample_df <- ParentModule$UncertainInput 
  uniqvNames <- unique(sample_df$varName)
  
  
  
  for (i in length(sample_df$data[1])){
    df <- sample_df |>
      select(varName, Scale, SubCompart)
    
    values <- first_values <- map(sample_df$data, ~ .x$value[i])
    
    df <- df |>
      mutate(Waarde = values)
    
    TheCore$mutateVars(df)
    
    #update core and solve
    TheCore$UpdateDirty(uniqvNames)
    ParentModule$PrepKaasM()
    ParentModule$PrepemisV()
  
    SB.K = ParentModule$SB.k
    vEmis = ParentModule$emissions
  
    solve(SB.K, -vEmis, tol = tol)
  }
  
  
  
  
  
}