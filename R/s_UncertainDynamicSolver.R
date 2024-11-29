#' @title UncertainDynamicSolver
#' @name UncertainDynamicSolver
#' @description solver to analyze sensitivity of parameters
#' @param ParentModule SBCore
#' @param tmax time [s] for the simulation period
#' @param sample_df nested tibble containing the sample values for every run and variable
#' @param nTIMES number of timesteps
#' @return Nested list containing the input variables, input emissions, output masses and states
#' @export
UncertainDynamicSolver = function(ParentModule, sample_df, nTIMES = 100,
                                  rtol_ode=1e-30, atol_ode = 1e-3) { 

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
  
  uniqvNames <- unique(sample_df$varName)
  
  # Get the emissions
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
  
  # Get the states
  StateAbbr <- rownames(ParentModule$SB.k)
  states <- ParentModule$myCore$states$asDataFrame |> 
    filter(Abbr %in% StateAbbr)
  
  # Create function to make approx functions from data (input is a df with the columns Abbr, Timed and Emis)
  makeApprox <- function(vEmissions){
    vEmis <- 
    vEmissions |> 
    group_by(Abbr) |> 
    summarise(n=n(),
              EmisFun = list(
                approxfun(
                  data.frame(Timed = c(0,Timed), 
                             Emis=c(0,Emis)),
                  rule = 1:2)
              )
    )
  funlist <- vEmis$EmisFun
  names(funlist) <- vEmis$Abbr
  return(funlist)
  } 
  
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
      Emis <- map_dfr(vEmissions$Emis, ~ tibble(Emis = .x$Mass_kg_s[i]))
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
    # browser()
    # Prepare the data to use mutateVar
    df <- sample_df |>
      unnest(data) |> filter(RUN == i) |> 
      select(varName, Scale, SubCompart, Species,value) |> 
      rename(Waarde = value)
    # values <- map(sample_df$data, ~ .x$value[i])
    # df <- df |>
    #   mutate(Waarde = values)
    ParentModule$myCore$mutateVars(df)

    # Update core
    ParentModule$myCore$UpdateDirty(uniqvNames)
    ParentModule$PrepKaasM()
    
    SB.K = ParentModule$SB.k
    
    SBNames = colnames(SB.K)
    SB.m0 <- rep(0, length(SBNames)) #TODO add this to input of solver.
    #SBtime <- seq(0,tmax,length.out = nTIMES)
    #SBtime <- seq(min(Emis_df$Timed), tmax, length.out = nTIMES)
    SBtime <- unique(Emis_df$Timed) # TODO this can be done differently
    
    # Define the solver function
    ODEapprox = function(t, m, parms) {
      with(as.list(c(parms, m)), {
        e <- c(rep(0, length(SBNames)))
        for (name in names(parms$emislist)) {
          e[grep(name, SBNames)] <- parms$emislist[[name]](t) 
        }
        dm <- with(parms, K %*% m + e) 
        return(list(dm, signal = e))
      })
    }
    
    # Solve the matrix
    sol <- deSolve::ode(
      y = as.numeric(SB.m0),
      times = SBtime,
      func = ODEapprox,
      parms = list(K = SB.K, SBNames=SBNames, emislist= funlist),
      rtol = rtol_ode, atol = atol_ode)
    
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