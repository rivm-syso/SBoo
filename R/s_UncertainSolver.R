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
  
  vEmissions = ParentModule$emissions
  
  states = ParentModule$myCore$states$asDataFrame
  
  SB.K = ParentModule$SB.k
  RowNames <- rownames(SB.K)
  states = states |>
    filter(Abbr %in% RowNames)
  
  for (i in 1:nrow(sample_df$data[[1]])){ 
    if(is.numeric(vEmissions$Emis)){
      Emis_df <- vEmissions
    } else if (inherits(vEmissions$Emis, "list")){
      Abbr <- vEmissions$Abbr
      Emis <- map_dfr(vEmissions$Emis, ~ tibble(Emis = .x$value[i]))
      Emis_df <- cbind(Abbr, Emis)
    }
    
    vEmis <- rep(0.0, length.out = length(RowNames))
    names(vEmis) <- states
    vEmis[match(Emis_df$Abbr, RowNames)] <- Emis_df$Emis
    names(vEmis) <- states
    
    df <- sample_df |>
      select(varName, Scale, SubCompart)
    
    values <- map(sample_df$data, ~ .x$value[i])
    
    df <- df |>
      mutate(Waarde = values)
    
    TheCore$mutateVars(df)
    
    #update core and solve
    TheCore$UpdateDirty(uniqvNames)
    
    ParentModule$PrepKaasM()

    SB.K = ParentModule$SB.k
    
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