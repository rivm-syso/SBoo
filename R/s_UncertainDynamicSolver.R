#' @title UncertainDynamicSolver
#' @name UncertainDynamicSolver
#' @description solver to analyze sensitivity of parameters
#' @param ParentModule SBCore
#' @param tmax time [s] for the simulation period
#' @param sample_df nested tibble containing the sample values for every run and variable
#' @param nTIMES number of timesteps
#' @return Nested list containing the input variables, input emissions, output masses and states
#' @export
UncertainDynamicSolver = function(ParentModule, tmax = 1e10, nTIMES = 100) { 
  
  # Get the uncertain input for the variables
  Input_Variables <- ParentModule$Input_Variables #including RUN + list

  uniqvNames <- unique(Input_Variables$varName)
  
  # Get the emissions
  vEmissions <- ParentModule$Emission$CleanEmissions
  
  for (i in 1:nrow(sample_df$data[[1]])){ 
    # First check: if vEmissions is a single emission data frame
    if(is.numeric(vEmissions$Emis)){ 
      
      funlist <- makeApprox(vEmissions)
      
    # Second check: if vEmissions is a single set of functions  
    } else if (inherits(vEmissions, "list")) {
      funlist <- vEmissions
      
    # Third check: if vEmissions is a nested tibble with emission data points  
    } else if (inherits(vEmissions$Emis, "list")){
      Abbr_Timed <- vEmissions |>
        select(Abbr, Timed)
      
      #Abbr <- vEmissions$Abbr
      Emis <- map_dfr(vEmissions$Emis, ~ tibble(Emis = .x$value[i]))
      Emis_df <- cbind(Abbr_Timed, Emis)
      
      funlist <- makeApprox(Emis_df)
      
    # Fourth check: if vEmissions is a nested tibble with sets of functions
    } else if ("Funlist" %in% colnames(vEmissions) && class(vEmissions$Funlist == "list")){
      
      subset_fun_tibble <- final_fun_tibble |>
        mutate(EmisFun = map(Funlist, ~ .x[[1]])) |>
        select(-Funlist)
      
      funlist <- subset_fun_tibble$EmisFun
      names(funlist) <- subset_fun_tibble$Abbr
    }
    
    # Prepare the data to use mutateVar
    df <- sample_df |>
      select(varName, Scale, SubCompart, Species)
    values <- map(sample_df$data, ~ .x$value[i])
    df <- df |>
      mutate(Waarde = values)
    ParentModule$myCore$mutateVars(df)
    
    # Update core
    ParentModule$myCore$UpdateDirty(uniqvNames)
    ParentModule$PrepKaasM()
    
    SB.K = ParentModule$SB.k
    
    SBNames = colnames(SB.K)
    SB.m0 <- rep(0, length(SBNames))
    #SBtime <- seq(0,tmax,length.out = nTIMES)
    #SBtime <- seq(min(Emis_df$Timed), tmax, length.out = nTIMES)
    SBtime <- unique(Emis_df$Timed)
    
    # Solve the matrix
    sol <- deSolve::ode(
      y = as.numeric(SB.m0),
      times = SBtime,
      func = ParentModule$ODEapprox,
      parms = list(K = SB.K, SBNames=SBNames, emislist= funlist),
      rtol = 1e-11, atol = 1e-3)
    
    # Change the colnames of the solved matrix
    colnames(sol)[1:length(SBNames)+1] <- SBNames
    colnames(sol)[grep("signal",colnames(sol))] <- paste("emis",SBNames,sep = "2")
    
    sol <- data.frame(sol)
    sol$RUN <- i
    
    if (!exists("solution")){
      solution <- sol
    } else {
      solution <- rbind(solution, sol)
    }
  }
  
  units <- World$fetchData("Units") |>
    select(VarName, Unit)
  
  sample_df <- left_join(sample_df, units, by = c("varName" = "VarName"))
  
  # if (inherits(vEmissions) == "data.frame") {
  #   vEmissions <- vEmissions |>
  #   mutate(Unit = "kg.s-1")
  # }
  # 
  # solution <- solution |>
  #   mutate(Unit = "kg")
  
  return(list(
    Input_Variables = sample_df,
    Input_Emission = vEmissions,
    DynamicMass = solution,
    States = states))
}