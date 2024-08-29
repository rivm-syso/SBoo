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
  sample_df <- ParentModule$UncertainInput 
  uniqvNames <- unique(sample_df$varName)
  
  solution <- sample_df |>
    select(varName, Scale, SubCompart) |>
    mutate(Waarde = 0) |>
    mutate(
      varName = ifelse(is.na(varName), "", varName),
      Scale = ifelse(is.na(Scale), "", Scale),
      SubCompart = ifelse(is.na(SubCompart), " ", SubCompart)
    ) %>%
    mutate(new_col_name = str_c(varName, Scale, SubCompart, sep = "_")) |> 
    select(new_col_name, Waarde) |>
    pivot_wider(names_from = new_col_name, values_from = Waarde) |>
    mutate(RUN = 0) |>
    mutate(Mass = 0)
  
  solution <- solution[-1,]
  
  l <- length(sample_df$data[1])
  
  vEmis = ParentModule$emissions
  
  for (i in 1:nrow(sample_df$data[[1]])){
    df <- sample_df |>
      select(varName, Scale, SubCompart)
    
    values <- first_values <- map(sample_df$data, ~ .x$value[i])
    
    df <- df |>
      mutate(Waarde = values)
    
    TheCore$mutateVars(df)
    
    #update core and solve
    TheCore$UpdateDirty(uniqvNames)
    #ParentModule$PrepKaasM()
    #preppedemis <- ParentModule$PrepemisV(emissions = ParentModule$emissions, solvername = "SBsteady") # Because we're solving for steady state, the emissions should be prepped accordingly
    
    SB.K = ParentModule$SB.k
    states = ParentModule$myCore$states$asDataFrame
    
    RowNames <- rownames(SB.K)
    states = states |>
      filter(Abbr %in% RowNames)
    
    sol <- solve(SB.K, -vEmis, tol=tol)
    
    sol <- tibble(sol) |>
      rename(EqMass = sol)
    
    sol <- cbind(states, sol)
    
    df <- df |> 
          mutate( 
            varName = ifelse(is.na(varName), "", varName),
            Scale = ifelse(is.na(Scale), "", Scale),
            SubCompart = ifelse(is.na(SubCompart), "", SubCompart)
          ) |>
          mutate(new_col_name = str_c(varName, Scale, SubCompart, sep = "_")) |> 
          select(new_col_name, Waarde) |>
          pivot_wider(names_from = new_col_name, values_from = Waarde) |>
          mutate(RUN = i) 
    
    sol_tibble <- tibble(Mass = list(sol))
    
    # Combine result_df and sol_tibble
    final_df <- bind_cols(df, sol_tibble)
      
    solution <- rbind(solution, final_df)
  }
  
  return(solution)
}