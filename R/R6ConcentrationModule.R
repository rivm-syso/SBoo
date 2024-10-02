#' @title ConcentrationModule
#' @description module for calculation of concentrations from masses
#' @param TheCore centre "World" object
#' @param input any other input that can be given, can be used to force the function to use a dataframe
#' @import R6
#' @export
ConcentrationModule <- 
  R6::R6Class("ConcentrationModule",
      public = list(
        initialize = function(TheCore, solvername){
          #browser()
          private$MyCore <- TheCore
          private$SolverName <- solvername
        
        },
        #'@description Leading function to obtain concentration from other world parameters
        #'@input all input is taken from the core
        GetConc = function(input){
          browser()
          if (private$SolverName == "SB1Solve" | private$SolverName == "SBsteady") {
           
            solution <- private$MyCore$Solution()
            states <- private$MyCore$states$asDataFrame
            volume <- private$MyCore$fetchData("Volume")
            fracw <- private$MyCore$fetchData("FRACw")
            fraca <- private$MyCore$fetchData("FRACa")
            rho <- private$MyCore$fetchData("rhoMatrix")
            
            longsolution <- solution |>
              left_join(states, by = "Abbr") |>
              left_join(private$MyCore$fetchData("Volume"), by = c("SubCompart", "Scale")) |>
              mutate(Concentration = EqMass / Volume)
            
            conctobecor <- longsolution  |>
              filter(SubCompart %in% c("marinesediment", "freshwatersediment", "lakesediment", 
                                       "agriculturalsoil", "naturalsoil", "othersoil"))
            conctobecor <- conctobecor |>
              left_join(fracw, by = c("SubCompart", "Scale")) |>
              left_join(fraca, by = c("SubCompart", "Scale")) |> 
              left_join(rho, by = c("SubCompart")) |> 
              mutate(rhowat = 998) #TODO: use rho from rhoMatrix
            
            conctobecor <- conctobecor |> 
              mutate(
                AdjustedConc = private$adjustconc(Concentration, FRACw, FRACa, rhoMatrix, rhowat)
              )
            
            longsolution <- longsolution |>
              left_join(
                conctobecor |> select(Scale, SubCompart, Species, AdjustedConc), 
                by = c("Scale", "SubCompart", "Species")
              ) |>
              mutate(
                Concentration = coalesce(AdjustedConc, Concentration)
              ) |>
              select(-AdjustedConc)
            
            subcompart <- c("agriculturalsoil", "naturalsoil", "othersoil", "freshwatersediment", "marinesediment",  "lakesediment", "air", "deepocean", "lake" , "river", "sea", "cloudwater")
            units <- c("g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw",
                       "g/m^3", "g/L", "g/L", "g/L", "g/L", "g/L", "g/L")
            
            # Combine into a named list
            subcompart_units <- setNames(units, subcompart)
            
            FinalConcentration <- longsolution |>
              select(Abbr, Scale, SubCompart, Concentration) |>
              mutate(Concentration = Concentration*1000) |>
              mutate(Unit = subcompart_units[SubCompart])
          
            return(FinalConcentration)
          }
          else if (private$SolverName == "SBsolve" | private$SolverName == "DynApproxSolve" | private$SolverName == "EventSolver") {
            #browser()
            solution <- private$MyCore$Solution()
            solution <- solution |>
              select(!starts_with("Emis"))
            rownames(solution) <- solution$time
            solution$time <- NULL 
            #transpose the data
            solution_t <- t(solution)
            
            #  Convert the transposed data back to a data frame
            solution_t <- as.data.frame(solution_t)
            
            #  Set the column names to be the transposed 'time' values
            colnames(solution_t) <- rownames(solution)
            solution_t$Abbr <- rownames(solution_t)
            #reorder
            solution_t <- solution_t[, c(ncol(solution_t), 1:(ncol(solution_t)-1))]
            
            
            states <- private$MyCore$states$asDataFrame
            volume <- private$MyCore$fetchData("Volume")
            
            longsolution <- solution_t |>
              left_join(states, by = "Abbr") |>
              left_join(private$MyCore$fetchData("Volume"), by = c("SubCompart", "Scale")) 
            
            #divide only the numeric (time) columns by Volume
            longsolution <- longsolution |> 
              mutate(across(
                .cols =matches(" ^\\d+$"),
                .fns = ~ . / Volume
              ))
            fracw <- private$MyCore$fetchData("FRACw")
            fraca <- private$MyCore$fetchData("FRACa")
            rho <- private$MyCore$fetchData("rhoMatrix")
            
            conctobecor <- longsolution  |>
              filter(SubCompart %in% c("marinesediment", "freshwatersediment", "lakesediment", 
                                       "agriculturalsoil", "naturalsoil", "othersoil"))
            conctobecor <- conctobecor |>
              left_join(fracw, by = c("SubCompart", "Scale")) |>
              left_join(fraca, by = c("SubCompart", "Scale")) |> 
              left_join(rho, by = c("SubCompart")) |> 
              mutate(rhowat = 998)
            
            conctobecor <- conctobecor |> 
              mutate(across(
                .cols = matches("^\\d+$"),  # Selects columns with numeric names
                .fns = ~ private$adjustconc(., FRACw, FRACa, rhoMatrix, rhowat)
              )) |>
              select(!c("Scale", "SubCompart", "Species", "Volume", "FRACw", "FRACa", "rhoMatrix", "rhowat"))
            common_cols <-setdiff(intersect(names(conctobecor), names(longsolution)), "Abbr")
            
            longsolution <- longsolution |>
              left_join(conctobecor, by = "Abbr", suffix = c("", ".new")) 
            
            longsolution <- longsolution |>
              mutate(across(all_of(common_cols), ~ coalesce(get(paste0(cur_column(), ".new")), .))) |>
              select(-ends_with(".new"))
            
            subcompart <- c("agriculturalsoil", "naturalsoil", "othersoil", "freshwatersediment", "marinesediment",  "lakesediment", "air", "deepocean", "lake" , "river", "sea", "cloudwater")
            units <- c("g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw",
                       "g/m^3", "g/L", "g/L", "g/L", "g/L", "g/L", "g/L")
            
            # Combine into a named list
            subcompart_units <- setNames(units, subcompart)
              
            FinalConcentration <- longsolution |>
              select(-Volume) |>
              mutate(across(where(is.numeric), ~ .x * 1000)) |>
              mutate(Unit = subcompart_units[SubCompart])
              
            
            return(FinalConcentration)
          } else if (private$SolverName == "UncertainSolver") {
            
            # TO DO: What if one of the variables below is uncertain?? 
            
            solution <- private$MyCore$Solution()$SteadyStateMass
            inputvars <- private$MyCore$Solution()$Input_Variables
            
            states <- private$MyCore$states$asDataFrame
            volume <- private$MyCore$fetchData("Volume")
            fracw <- private$MyCore$fetchData("FRACw")
            fraca <- private$MyCore$fetchData("FRACa")
            rho <- private$MyCore$fetchData("rhoMatrix")
            
            longsolution <- solution |>
              left_join(private$MyCore$fetchData("Volume"), by = c("SubCompart", "Scale")) |>
              mutate(Concentration = EqMass / Volume)
            
            conctobecor <- longsolution  |>
              filter(SubCompart %in% c("marinesediment", "freshwatersediment", "lakesediment", 
                                       "agriculturalsoil", "naturalsoil", "othersoil"))
            conctobecor <- conctobecor |>
              left_join(fracw, by = c("SubCompart", "Scale")) |>
              left_join(fraca, by = c("SubCompart", "Scale")) |> 
              left_join(rho, by = c("SubCompart")) |> 
              mutate(rhowat = 998) #TODO: use rho from rhoMatrix
            
            conctobecor <- conctobecor |> 
              mutate(
                AdjustedConc = private$adjustconc(Concentration, FRACw, FRACa, rhoMatrix, rhowat)
              )
            
            longsolution <- longsolution |>
              left_join(
                conctobecor |> select(Scale, SubCompart, Species, RUN, AdjustedConc), 
                by = c("Scale", "SubCompart", "Species", "RUN")
              ) |>
              mutate(
                Concentration = coalesce(AdjustedConc, Concentration)
              ) |>
              select(-AdjustedConc)
            
            subcompart <- c("agriculturalsoil", "naturalsoil", "othersoil", "freshwatersediment", "marinesediment",  "lakesediment", "air", "deepocean", "lake" , "river", "sea", "cloudwater")
            units <- c("g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw",
                       "g/m^3", "g/L", "g/L", "g/L", "g/L", "g/L", "g/L")
            
            # Combine into a named list
            subcompart_units <- setNames(units, subcompart)
            
            FinalConcentration <- longsolution |>
              select(Abbr, Scale, SubCompart, RUN, Concentration) |>
              mutate(Concentration = Concentration*1000) |>
              mutate(Unit = subcompart_units[SubCompart])
            
            return(FinalConcentration)
          } else if (private$SolverName == "UncertainDynamicSolver") {
            
            # TO DO: Adjust this code for the df that includes RUN, most is copied from normal dynamic solver concentrations above. 
            
            solution <- private$MyCore$Solution()$DynamicMass |>
              select(!starts_with("emis"))
            rownames(solution) <- solution$time
            solution$time <- NULL 
            
            #transpose the data
            solution_t <- t(solution)
            
            #  Convert the transposed data back to a data frame
            solution_t <- as.data.frame(solution_t)
            
            #  Set the column names to be the transposed 'time' values
            colnames(solution_t) <- rownames(solution)
            solution_t$Abbr <- rownames(solution_t)
            #reorder
            solution_t <- solution_t[, c(ncol(solution_t), 1:(ncol(solution_t)-1))]
            
            inputvars <- private$MyCore$Solution()$Input_Variables
            
            # TO DO: What if one of the variables below is uncertain?? 
            states <- private$MyCore$states$asDataFrame
            volume <- private$MyCore$fetchData("Volume")
            fracw <- private$MyCore$fetchData("FRACw")
            fraca <- private$MyCore$fetchData("FRACa")
            rho <- private$MyCore$fetchData("rhoMatrix")
  
            states <- private$MyCore$states$asDataFrame

            longsolution <- solution_t |>
              left_join(states, by = "Abbr") |>
              left_join(private$MyCore$fetchData("Volume"), by = c("SubCompart", "Scale")) 
            
            #divide only the numeric (time) columns by Volume
            longsolution <- longsolution |> 
              mutate(across(
                .cols =matches(" ^\\d+$"),
                .fns = ~ . / Volume
              ))

            conctobecor <- longsolution  |>
              filter(SubCompart %in% c("marinesediment", "freshwatersediment", "lakesediment", 
                                       "agriculturalsoil", "naturalsoil", "othersoil"))
            conctobecor <- conctobecor |>
              left_join(fracw, by = c("SubCompart", "Scale")) |>
              left_join(fraca, by = c("SubCompart", "Scale")) |> 
              left_join(rho, by = c("SubCompart")) |> 
              mutate(rhowat = 998)
            
            conctobecor <- conctobecor |> 
              mutate(across(
                .cols = matches("^\\d+$"),  # Selects columns with numeric names
                .fns = ~ private$adjustconc(., FRACw, FRACa, rhoMatrix, rhowat)
              )) |>
              select(!c("Scale", "SubCompart", "Species", "Volume", "FRACw", "FRACa", "rhoMatrix", "rhowat"))
            common_cols <-setdiff(intersect(names(conctobecor), names(longsolution)), "Abbr")
            
            longsolution <- longsolution |>
              left_join(conctobecor, by = "Abbr", suffix = c("", ".new")) 
            
            longsolution <- longsolution |>
              mutate(across(all_of(common_cols), ~ coalesce(get(paste0(cur_column(), ".new")), .))) |>
              select(-ends_with(".new"))
            
            subcompart <- c("agriculturalsoil", "naturalsoil", "othersoil", "freshwatersediment", "marinesediment",  "lakesediment", "air", "deepocean", "lake" , "river", "sea", "cloudwater")
            units <- c("g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw",
                       "g/m^3", "g/L", "g/L", "g/L", "g/L", "g/L", "g/L")
            
            # Combine into a named list
            subcompart_units <- setNames(units, subcompart)
            
            FinalConcentration <- longsolution |>
              select(-Volume) |>
              mutate(across(where(is.numeric), ~ .x * 1000)) |>
              mutate(Unit = subcompart_units[SubCompart])
            
          } else { 
            stop("Can not calculate concentration for this solver")
          }
        }
      ),
        # concentrationcalc = function() 
      private = list(
        MyCore = NULL, 
        SolverName = NULL,
        #'@name Adjustconc
        #'@description helper function to adjust concentration for soil and sediment. Returns kg/kg dry weight.
        #'@param CompConc Concentration in each compartment [mass/volume]
        #'@param FRACw the fraction of water in the matrix
        #'@param FRACa the fraction of air in the matrix 
        #'@param RHOsolid the density of the solid phase of the material [kg m-3]
        #'@param RhoWater_value the density of water [kg m-3]
        adjustconc = function(CompConc, Fracw, Fraca, RHOsolid, RhoWater_value){
            CompConc  / (Fracw * RhoWater_value + (1 - Fracw - Fraca) * RHOsolid)
          }
      )
  )

