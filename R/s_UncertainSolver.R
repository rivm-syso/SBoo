#' @title UncertainSolver
#' @name UncertainSolver
#' @description solver to include uncertainty and variability of input variables
#' @param ParentModule SBCore
#' @param sample_df nested tibble containing the sample values for every run and variable
#' @param n samplesize 
#' @return States (i) (=mass)
#' @export
UncertainSolver = function(ParentModule, tol=1e-10, sample_df) { 
  
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
  
  # Get the emissions and states
  vEmissions = ParentModule$emissions
  
  if(!is.numeric(vEmissions$Emis)){
    if(all(map_lgl(vEmissions$Emis, ~ "RUN" %in% names(.x))) == FALSE){
      warning("adding RUN number to emission data")
      vEmissions <- vEmissions |> 
        mutate(nRUNs = map_int(Emis, nrow)) |> 
        mutate(
          Emis = map(Emis, ~ .x |> 
                       mutate(RUN = 1:unique(nRUNs)))
        ) |> select(-nRUNs)
    }
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

    VariableInputRun <- sample_df |> 
      mutate(Waarde =  map_vec(sample_df$data, ~ .x$value[i])) |> 
      select(-data)

    ParentModule$myCore$mutateVars(VariableInputRun)
    
    #update core and solve
    ParentModule$myCore$UpdateDirty(unique(VariableInputRun$varName))
    ParentModule$PrepKaasM()
    
    if(det(ParentModule$SB.k == 0)){
      sol <- ginv(ParentModule$SB.k) %*% (-vEmis)
      print("K matrix is singular, generalized inverse was used to solve the matrix")
    } else {
      sol <- solve(ParentModule$SB.k, # K matrix of first order rate constants
                 -vEmis, # Emission vector
                 tol=tol) # solve tolerance
    }
    
    sol <- tibble(EqMass = sol, 
                  Abbr = names(sol),
                  RUN = i) |> 
      full_join(states)
    
    if(!exists("solution")) solution <- data.frame(NULL) # create on first loop
    solution <- rbind(solution, sol)
  }
  
  units <- World$fetchData("Units") |>
    select(VarName, Unit)
  
  sample_df <- left_join(sample_df, units, by = c("varName" = "VarName"))
  
  vEmissions <- vEmissions |>
    mutate(Unit = "kg.s-1")
  
  solution <- solution |>
    mutate(Unit = "kg")
  
  return(list(
    Input_Variables = sample_df,
    Input_Emission = vEmissions,
    SteadyStateMass = solution))
}
