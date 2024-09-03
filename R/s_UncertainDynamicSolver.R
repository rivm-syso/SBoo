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
  nt <- nTIMES
  
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
    
    # Prepare the emissions like in the emission module
    
    # Emission preparation for when a df with timed column is given
    vEmis <- 
      Emis_df |> 
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
    
    # Prepare the data to use mutateVar
    df <- sample_df |>
      select(varName, Scale, SubCompart)
    values <- map(sample_df$data, ~ .x$value[i])
    df <- df |>
      mutate(Waarde = values)
    ParentModule$myCore$mutateVars(df)
    
    #Update core
    ParentModule$myCore$UpdateDirty(uniqvNames)
    ParentModule$PrepKaasM()
    
    SB.K = ParentModule$SB.k
    
    SBNames = colnames(SB.K)
    SB.m0 <- rep(0, length(SBNames))
    SBtime <- seq(0,tmax,length.out = nTIMES)
    
    ODEapprox = function(t, m, parms) {
      #browser()
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
    
    colnames(sol)[1:length(SBNames)+1] <- SBNames
    colnames(sol)[grep("signal",colnames(sol))] <- paste("emis",SBNames,sep = "2")
    #new_colnames <- colnames(sol)
    
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
    
    #sol_tibble <- tibble(sol) |>
    #  nest(Mass = everything())
    sol <- tibble(sol) 
    
    sol_tibble <- tibble(Mass = list(sol))
    
    # Combine result_df and sol_tibble
    final_df <- bind_cols(df, sol_tibble)
    
    solution <- rbind(solution, final_df)
  }
  
  return(solution)
}