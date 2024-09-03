#' @title UncertainSolver
#' @name UncertainSolver
#' @description solver to analyze sensitivity of parameters
#' @param ParentModule SBCore
#' @param vnamesDistSD  dataframe with columns vnames, distNames (see lhs package for possible distributions), secondPar
#' @param n samplesize 
#' @return States (i) (=mass)
#' @export
UncertainSolver = function(ParentModule, tol=1e-30) { 
  
  # Get the uncertain input for the variables
  sample_df <- ParentModule$UncertainInput 
  if(all(map_lgl(sample_df$data, ~ "RUN" %in% names(.x))) == FALSE){
    warning("adding RUN number to variable data")
    sample_df <- sample_df |> 
      mutate(nRUNs = map_int(data, nrow)) |> 
      mutate(
        data = map(data, ~ .x |> 
                     mutate(RUN = 1:unique(nRUNs)))
      ) |> select(-nRUNs)
    
  }
  
  # Create an empty df to store the final solution in after the for loop
  # solution <- sample_df |>
  #   select(varName, Scale, SubCompart) |>
  #   mutate(Waarde = 0) |>
  #   mutate(
  #     varName = ifelse(is.na(varName), "", varName),
  #     Scale = ifelse(is.na(Scale), "", Scale),
  #     SubCompart = ifelse(is.na(SubCompart), " ", SubCompart)
  #   ) %>%
  #   mutate(new_col_name = str_c(varName, Scale, SubCompart, sep = "_")) |> 
  #   select(new_col_name, Waarde) |>
  #   pivot_wider(names_from = new_col_name, values_from = Waarde) |>
  #   mutate(RUN = 0) |>
  #   mutate(Mass = 0)
  # solution <- solution[-1,]
  
  # Get the emissions and states
  vEmissions = ParentModule$emissions
  if(all(map_lgl(vEmissions$Emis, ~ "RUN" %in% names(.x))) == FALSE){
    warning("adding RUN number to emission data")
    vEmissions <- vEmissions |> 
      mutate(nRUNs = map_int(Emis, nrow)) |> 
      mutate(
        Emis = map(Emis, ~ .x |> 
                     mutate(RUN = 1:unique(nRUNs)))
      ) |> select(-nRUNs)
    
  }
  
  StateAbbr <- rownames(ParentModule$SB.k)
  states <- ParentModule$myCore$states$asDataFrame |> 
    filter(Abbr %in% StateAbbr)
  
  #TODO make a for-each loop with ncores as variable.
  
  for (i in 1:nrow(sample_df$data[[1]])){ 
    if(is.numeric(vEmissions$Emis)){
      Emis_df <- vEmissions
    } else if (inherits(vEmissions$Emis, "list")){
      Abbr <- vEmissions$Abbr
      Emis <- map_dfr(vEmissions$Emis, ~ tibble(Emis = .x$value[i]))
      Emis_df <- cbind(Abbr, Emis)
    } else stop("no vEmissions found")
    
    vEmis <- rep(0.0, length.out = length(StateAbbr))
    names(vEmis) <- states$Abbr
    vEmis[match(Emis_df$Abbr, StateAbbr)] <- Emis_df$Emis
    # names(vEmis) <- states
    
    VariableInputRun <- sample_df |> 
      mutate(Waarde =  map_vec(sample_df$data, ~ .x$value[i])) |> 
      select(-data)
    #   select(varName, Scale, SubCompart)
    # 
    # values <- map_vec(sample_df$data, ~ .x$value[i])
    # VariableInputRun <- VariableInputRun |>
    #   mutate(Waarde = values)
    
    ParentModule$myCore$mutateVars(VariableInputRun)
    
    #update core and solve
    ParentModule$myCore$UpdateDirty(unique(VariableInputRun$varName))
    ParentModule$PrepKaasM()
    
    sol <- solve(ParentModule$SB.k, # K matrix of first order rate constants
                 -vEmis, # Emission vector
                 tol=tol) # solve tolerance
    
    sol <- tibble(EqMass = sol, 
                  Abbr = names(sol),
                  RUN = i) |> 
      full_join(states)
    
    # df <- 
    #   df |> 
    #   mutate( 
    #     varName = ifelse(is.na(varName), "", varName),
    #     Scale = ifelse(is.na(Scale), "", Scale),
    #     SubCompart = ifelse(is.na(SubCompart), "", SubCompart)
    #   ) |>
    #   mutate(new_col_name = str_c(varName, Scale, SubCompart, sep = "_")) |> 
    #   select(new_col_name, Waarde) |>
    #   pivot_wider(names_from = new_col_name, values_from = Waarde) |>
    #   mutate(RUN = i) 
    # 
    # sol_tibble <- tibble(Mass = list(sol))
    # 
    # # Combine result_df and sol_tibble
    # final_df <- bind_cols(df, sol_tibble)
    if(!exists("solution")) solution <- data.frame(NULL) # create on first loop
    solution <- rbind(solution, sol)
  }
  
  return(list(
    Input_Variables = sample_df,
    Input_Emission = vEmissions,
    SteadyStateMass = solution))
}
