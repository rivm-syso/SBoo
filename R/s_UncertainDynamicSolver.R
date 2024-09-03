#' @title UncertainDynamicSolver
#' @name UncertainDynamicSolver
#' @description solver to analyze sensitivity of parameters
#' @param ParentModule SBCore
#' @param vnamesDistSD  dataframe with columns vnames, distNames (see lhs package for possible distributions), secondPar
#' @param n samplesize 
#' @return States (i) (=mass)
#' @export
UncertainDynamicSolver = function(ParentModule, tmax = 1e10, sample_df, nTIMES = 100) { 
  
  sample_df <- ParentModule$UncertainInput 
  uniqvNames <- unique(sample_df$varName)
 
  # Create an empty tibble to bind the solutions to in the for loop
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
  
  # Filter the states to only contain the ones that are in the colnames of SB.K
  SB.K = ParentModule$SB.k
  RowNames <- rownames(SB.K)
  states = states |>
    filter(Abbr %in% RowNames)
  
  for (i in 1:nrow(sample_df$data[[1]])){ 
    # First check: if vEmissions is a single emission data frame
    if(is.numeric(vEmissions$Emis)){                
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
      
    # Second check: if vEmissions is a single set of functions  
    } else if (class(vEmissions) == "list") {
      vEmis <- vEmissions
      
    # Third check: if vEmissions is a nested tibble with emission data points  
    } else if (inherits(vEmissions$Emis, "list")){
      Abbr <- vEmissions$Abbr
      Emis <- map_dfr(vEmissions$Emis, ~ tibble(Emis = .x$value[i]))
      Emis_df <- cbind(Abbr, Emis)
      
    # Fourth check: if vEmissions is a nested tibble with sets of functions
    } #else if {
      
    #}
    
    # Prepare the data to use mutateVar
    df <- sample_df |>
      select(varName, Scale, SubCompart)
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
    SBtime <- seq(0,tmax,length.out = nTIMES)
    
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
      times = SBtime ,
      func = ODEapprox,
      parms = list(K = SB.K, SBNames=SBNames, emislist= funlist),
      rtol = 1e-11, atol = 1e-3)
    
    # Change the colnames of the solved matrix
    colnames(sol)[1:length(SBNames)+1] <- SBNames
    colnames(sol)[grep("signal",colnames(sol))] <- paste("emis",SBNames,sep = "2")
    
    # Make a tibble in the same format as the 'Solution' tibble
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
    
    # Make a tibble from the solved matrix and nest it in a column names Mass
    sol <- tibble(sol) 
    sol_tibble <- tibble(Mass = list(sol))
    
    # Combine result_df and sol_tibble
    final_df <- bind_cols(df, sol_tibble)
    
    # Bind the new solution to the previous solutions in the solution matrix
    solution <- bind_rows(solution, final_df)
  }
  
  return(solution)
}