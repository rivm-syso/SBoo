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
          #browser()
          if (private$SolverName == "SB1Solve" | private$SolverName == "SBsteady") {
           
            solution <- private$MyCore$Solution()
            states <- private$MyCore$states$asDataFrame
            
            # Add air and cloudwater volume together
            volume <- private$MyCore$fetchData("Volume") |>
              mutate(SubCompart = ifelse(SubCompart == "cloudwater", "air", SubCompart)) |>
              group_by(Scale, SubCompart) |>
              summarise(Volume = sum(Volume)) |>
              ungroup()
            fracw <- private$MyCore$fetchData("FRACw")
            fraca <- private$MyCore$fetchData("FRACa")
            rho <- private$MyCore$fetchData("rhoMatrix")
            
            # sum mass in air and cloudwater compartments and calculate the concentration
            longsolution <- solution |>
              left_join(states, by = "Abbr") |>
              mutate(SubCompart = ifelse(SubCompart == "cloudwater", "air", SubCompart)) |>
              mutate(Abbr = str_replace_all(Abbr, "cw", "a")) |>
              group_by(Scale, SubCompart, Species, Abbr) |>
              summarise(EqMass = sum(EqMass)) |>
              ungroup() |>
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
            
            subcompart <- c("agriculturalsoil", "naturalsoil", "othersoil", "freshwatersediment", "marinesediment",  "lakesediment", "air", "deepocean", "lake" , "river", "sea")
            units <- c("g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw",
                       "g/m^3", "g/L", "g/L", "g/L", "g/L", "g/L")
            
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
            volume <- private$MyCore$fetchData("Volume") |>
              mutate(SubCompart = ifelse(SubCompart == "cloudwater", "air", SubCompart)) |>
              group_by(Scale, SubCompart) |>
              summarise(Volume = sum(Volume)) |>
              ungroup()
            
            longsolution <- solution_t |>
              left_join(states, by = "Abbr") 
            
            longsolution <- longsolution |>
              mutate(SubCompart = ifelse(SubCompart == "cloudwater", "air", SubCompart)) |>
              mutate(Abbr = str_replace_all(Abbr, "cw", "a")) |>
              group_by(Scale, SubCompart, Species, Abbr) |>
              summarise(across(where(~ is.numeric(.x)), sum, na.rm = TRUE)) |>
              ungroup() |>
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
            
            subcompart <- c("agriculturalsoil", "naturalsoil", "othersoil", "freshwatersediment", "marinesediment",  "lakesediment", "air", "deepocean", "lake" , "river", "sea")
            units <- c("g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw", "g/kg dw",
                       "g/m^3", "g/L", "g/L", "g/L", "g/L", "g/L")
            
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
            volume <- private$MyCore$fetchData("Volume") |>
              mutate(SubCompart = ifelse(SubCompart == "cloudwater", "air", SubCompart)) |>
              group_by(Scale, SubCompart) |>
              summarise(Volume = sum(Volume)) |>
              ungroup()
            fracw <- private$MyCore$fetchData("FRACw")
            fraca <- private$MyCore$fetchData("FRACa")
            rho <- private$MyCore$fetchData("rhoMatrix")
            
            longsolution <- solution |>
              mutate(SubCompart = ifelse(SubCompart == "cloudwater", "air", SubCompart)) |>
              mutate(Abbr = str_replace_all(Abbr, "cw", "a")) |>
              group_by(RUN, Scale, SubCompart, Species, Unit, Abbr) |>
              summarise(EqMass = sum(EqMass)) |>
              ungroup() |>
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
            
            solution <- private$MyCore$Solution()$DynamicMass |>
              select(!starts_with("emis"))
            
            # sum the cw values with the a values in the same scale/species
            cw_cols <- grep("^cw", colnames(solution), value = TRUE)
            
            # Loop through each "cw" column
            for (cw_col in cw_cols) {
              
              # Extract the suffix after "cw"
              suffix <- sub("^cw", "", cw_col)
              
              # Construct the corresponding "a" column name
              a_col <- paste0("a", suffix)
              
              # Check if the "a" column exists
              if (a_col %in% colnames(solution)) {
                # Add the "cw" column values to the "a" column
                solution[[a_col]] <- solution[[a_col]] + solution[[cw_col]]
              }
            }
            
            # TO DO: What if one of the variables below is uncertain?? 
            states <- private$MyCore$states$asDataFrame
            volume <- private$MyCore$fetchData("Volume") |>
              mutate(SubCompart = ifelse(SubCompart == "cloudwater", "air", SubCompart)) |>
              group_by(Scale, SubCompart) |>
              summarise(Volume = sum(Volume)) |>
              ungroup()
            fracw <- private$MyCore$fetchData("FRACw")
            fraca <- private$MyCore$fetchData("FRACa")
            rho <- private$MyCore$fetchData("rhoMatrix")

            compnames <- colnames(solution)
            compnames <- compnames[compnames %in% states$Abbr]
            
            varvalues <- states |>
              left_join(volume, by = c("Scale", "SubCompart")) |>
              left_join(fracw, by = c("Scale", "SubCompart")) |>
              left_join(fraca, by = c("Scale", "SubCompart")) |>
              left_join(rho, by = "SubCompart") |>
              select(!c(Scale, SubCompart, Species)) |>
              filter(Abbr %in% compnames) |>
              mutate(rhowat  = 998)           # Add the density of water, needed later
              
            varvalues_t <- t(varvalues)
            varvalues_t <- data.frame(varvalues_t)
            
            colnames(varvalues_t) <- varvalues$Abbr
            
            varvalues_t <- varvalues_t[-1, ]
            varvalues_t[] <- lapply(varvalues_t, as.numeric)
            
            nruns_sol <- nrow(solution)
            nrow_var <- nrow(varvalues_t)
            
            solution <- full_join(solution, varvalues_t, by = colnames(varvalues_t))
            
            rownames(solution)[(nruns_sol+1):(nruns_sol+nrow_var)] <- rownames(varvalues_t)
            
            inputvars <- private$MyCore$Solution()$Input_Variables
            
            # Make a list of compnames that should be corrected later (all sediment and soil compartments)
            comptobecor <- grep("^s", compnames, value = TRUE)
            
            comptokeep <- grep("^s", compnames, value = TRUE, invert = TRUE)
            
            concentration_df <- solution |>
              mutate(time = as.character(time)) |>
              mutate(RUN = as.character(RUN))

            numeric_cols <- concentration_df |>
              select(all_of(compnames))
            
            other_cols <- concentration_df |>
              select(-all_of(compnames))
            
            volume_values <- as.numeric(numeric_cols["Volume", ])
            
            concentrations <- apply(numeric_cols[1:nruns_sol, ], 2, function(column) {
              (column / volume_values) *1000
            }) 
            
            #concentration_df[1:nruns_sol, ] <- concentrations
            concentrations <- data.frame(concentrations) 
            
            concentrations <- full_join(concentrations, varvalues_t, by = colnames(varvalues_t))
            concentrations <- cbind(concentrations, other_cols)
            
            conctobecor <- concentrations |>
              select(all_of(comptobecor))
            
            # Save the compartments that have the right concentrations already
            conctokeep <- concentrations |>
              select(time, RUN, all_of(comptokeep)) 
            
            # Calculate the concentrations for the compartments that need correcting
            FRACw <- conctobecor["FRACw", ]
            FRACa <- conctobecor["FRACa", ]
            rhoMatrix <- conctobecor["rhoMatrix", ]
            rhowat <- conctobecor["rhowat", ]
            
            # Apply the adjustment using mapply to pass multiple arguments to the function
            conctobecor_adj <- as.data.frame(mapply(function(column, FRACw, FRACa, rhoMatrix, rhowat) {
              private$adjustconc(column, FRACw, FRACa, rhoMatrix, rhowat)
            }, conctobecor[1:nruns_sol, ], FRACw, FRACa, rhoMatrix, rhowat))
            
            Units <- states |>
              mutate(Unit = case_when(
                startsWith(Abbr, "s") ~ "g/kg dw",
                startsWith(Abbr, "w") ~ "g/L",
                startsWith(Abbr, "a") ~ "g/m^3",
                TRUE ~ NA
              )) |>
              select(Abbr, Unit) |>
              filter(Abbr %in% compnames)
              
            Units <- t(Units)
            Units <- data.frame(Units, stringsAsFactors = FALSE)
            
            colnames(Units) <- Units[1, ]
            Units <- Units[-1, ]
            colnames(Units) <- make.names(colnames(Units), unique = TRUE)

            conctokeep <- conctokeep |>
              slice(1:(n()-nrow_var))
            
            final_concentrations <- cbind(conctokeep, conctobecor_adj)
            
            return(list(
              Concentrations = final_concentrations,
              Units = Units))
              
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

