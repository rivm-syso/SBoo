#' @title SolverModule
#' @description
#' wrapper for your solver function and collection of general methods needed for most solvers,
#' like preparing K matrix and matching emission vector and
#' the general ODE for SimpleBox
#' @import R6
#' @export
SolverModule <-
  R6::R6Class(
    "SolverModule",
    inherit = CalcGraphModule,

    private = list(
      NeedVars = function() {
        #overrule CalcGraphModule
        NULL
      },
      Masses = NULL,
      Concentration = NULL,
      UsedEmissions = NULL,
      AllVars = NULL,
      SolveStates = NULL,
      emissionModule = NULL,
      SB.K = NULL,
      MatrixMassesInRows = NULL,
      lSBtime.tvars = NULL,
      lvnamesDistSD = NULL,
      input_variables = NULL,
      ConcParams = NULL,
      Mass2ConcFun = NULL,
      run_df = NULL,
      LHSruns = matrix() # LFS matrix for the uncertainty variables
      ),
    
    public = list(
      initialize = function(TheCore, exeFunction, ...) {
        # Call the parent class's initialize method; inits also private$MyName
        super$initialize(TheCore, exeFunction, ...)
        
        # a SB variable outside of the ModuleList and nodelist of the core
        private$Mass2ConcFun <-
          VariableModule$new(TheCore, "Mass2ConcDivider", IsVectorised = T, AggrBy = NULL, AggrFun = NULL)

      },
      #' @description Solve matrix / emissions
      #' @param needdebug if T, the solver defining function execute in debugmode
      #' @param emissions
      #' @param var_box_df = data.frame(), 
      #' @param var_invFun = list(), 
      #' @param emis_invFun = list(), 
      #' @param nRUNs = NULL
      execute = function(needdebug = F,
                         emissions = NULL,
                         var_box_df = data.frame(), var_invFun = list(), nRUNs = NULL, 
                         ParallelPreparation=F, LHSmatrix = NULL, 
                         ...) {

        lhsRUNS <- 0 #unless overwritten:
        nVars <- length(var_invFun)
        
        if(is.null(nRUNs) || is.na(nRUNs)) {
          nRUNs <- 1
        }
        
        # If ParallelPreparation is TRUE, the LHS samples and World need to be saved for later use. 
        if(ParallelPreparation == T){
          
          # Prep LHS samples
          if (nRUNs < 2){
            stop("nRUNs should be greater than 1")
          } else if(length(var_invFun) == 0){
            stop("Please provide inverse functions")
          } else if(nrow(var_box_df) == 0){
            stop("Please provide variable dataframe")
          } else if (nRUNs != 0 && nRUNs != 1) {
            if (nVars > 0){
              # create the RUNs for variables and possibly emissions
              # NB PrepLHS sets private$input_variables to var_box_df
              lhsRUNS <- self$PrepLHS(var_box_df, var_invFun, emis_invFun = NULL, nRUNs)

              scaled_samples <- self$ScaleLHS(lhsRUNS, var_invFun)
              
              private$LHSruns <- t(scaled_samples)
              # Mind the transpose, for easy substituting the samples
              idnames <- c("varName", names(var_box_df)[names(var_box_df) %in% The3D])
              rownames(private$LHSruns) <- do.call(paste, as.list(var_box_df[,idnames]))
              
              saveRDS(private$LHSruns, "data/scaledLHSsamples.RData")
            }
          }
            # Export World
            saveRDS(World, file = "data/World.RData")
        
        # If ParallelPreparation is FALSE, regular solver use.
        } else {
          if (is.null(private$SB.K)) {
            stop("run PrepKaasM() first") #solveStates should be known, set by PrepKaasM()
          }
          
          MoreParams <- list(...)
          if (length(MoreParams) != 0) {
            nTIMES <- MoreParams$nTIMES
            tmax <- MoreParams$tmax
          } else {
            nTIMES <- NULL
            tmax <- NULL
          }
          
          if (nRUNs > 1) {
            # Create a matrix with original_runs and solver_runs
            used_runs <- 1:nRUNs
            solver_runs <- unique(emissions$RUN)
            
            run_matrix <- cbind(used_runs, solver_runs)
            private$run_df <- as.data.frame(run_matrix)
            
            emissions$RUN <- private$run_df$used_runs[match(emissions$RUN, private$run_df$solver_runs)]
            
          } else {
            # Create a matrix with a single row
            private$run_df <- data.frame(used_runs = 1,
                                         solver_runs = 1)
            
          }
          
          # Check if all elements in 'emissions' are functions
          if (all(sapply(emissions, is.function))) {
            emisFuns <- length(emissions)  # Number of functions in 'emissions'
          } else {
            emisFuns <- 0  # Not all elements are functions
          }
          
          # are the functions inv.lhs or in time? # or use nTIMES as criterion?
          # Determine the value of 'argName' based on 'emisFuns'
          if (emisFuns == 0) {
            argName <- ""
          } else {
            argName <- unique(sapply(emissions, formalArgs))
          }
          
          # Ensure that 'argName' has exactly one element
          stopifnot(length(argName) == 1)
          
          # This means the functions for emissions are approxfuns, this imput is not implemented currently
          if (argName %in% c("t", "v")) {
            
            stop("Approx funs as emissions are not implemented yet")
            
          } else if(argName == ""){
            
            self$PrepemisV(emissions)
            nEmisComps <- 0 # prevent emissions from being overwritten - no LHS samples will be taken for emissions
            
          } else { # This means the functions for emissions are distribution functions: for steady solver only.
            # we need to prep lhs; possibly combined with uncertainty in variables. If so, we need to know how many compartments have a distribution function
            nEmisComps <- length(emissions)
          }
          
          # Check if LHS samples still need to be prepared
          if(is.null(LHSmatrix)){
            # If nRUNs is specified by the user, pull samples from emission/variable distributions
            if (nRUNs != 0 && nRUNs != 1) {
              if (nVars + nEmisComps > 0){
                # create the RUNs for variables and possibly emissions
                # NB PrepLHS sets private$input_variables to var_box_df
                if (nEmisComps > 0) {
                  lhsRUNS <- self$PrepLHS(var_box_df, var_invFun, emis_invFun = emissions, nRUNs)
                } else {
                  lhsRUNS <- self$PrepLHS(var_box_df, var_invFun, emis_invFun = NULL, nRUNs)
                }
              }
            } else { 
              self$PrepemisV(emissions)
            }
            
            # split the lhs results over vars and emissions, apply to emissions
            if (nEmisComps > 1) {
              emissRuns <- lhsRUNS[,1:nEmisComps]
              colnames(emissRuns) <- names(emissions)
              self$PrepemisV(emissRuns)
            }
            if (nVars > 0) { # MIND YOU: the private lhsruns are for variables only!
              varLHS <- lhsRUNS[,(nEmisComps+1):ncol(lhsRUNS)]
              scaled_samples <- self$ScaleLHS(varLHS, var_invFun)
              inputvars <- private$input_variables
              private$LHSruns <- t(scaled_samples)
              # Mind the transpose, for easy substituting the samples
              idnames <- c("varName", names(var_box_df)[names(var_box_df) %in% The3D])
              rownames(private$LHSruns) <- do.call(paste, as.list(var_box_df[,idnames]))
            } else {
              nVars = 1
            }
          # If LHS samples are already prepared, assign LHSmatrix to private$LHSruns  
          } else {

            private$LHSruns <- LHSmatrix
            
            # Get the number of runs
            nRUNs <- ncol(private$LHSruns)
            
            # Get the number of variables
            nVars <- dim(private$LHSruns)[1] 
            
            # Set the input_variables dataframe from the matrix
            df <- as.data.frame(LHSmatrix)
            df <- cbind(do.call(rbind, strsplit(rownames(df), " ")), df)
            colnames(df)[1:4] <- c("varName", "Scale", "SubCompart", "Species")
            rownames(df) <- 1:nrow(df)
            df <- df[1:4]
            df[df == "NA"] <- NA
            
            private$input_variables <- df
            
            # Check if the length of the emission df is the same as nRUNs
            nRUNsEmis <- length(unique(emissions$RUN))
            if(!nRUNsEmis == nRUNs){
              stop("Different number of runs in LHS samples and emissions")
            }
          }
          
          if(is.null(nTIMES)){
            nTIMES <- 1
          } 
          
          if(is.null(tmax)){
            tmax = 0
          }
          
          # the resulting array is (allocated once)
          private$Masses <- array(dim = c(nTIMES,self$solveStates$nStates, nRUNs),
                                  dimnames = list(
                                    time = seq(0, tmax, length.out = nTIMES),
                                    self$solveStates$asDataFrame$Abbr,
                                    RUNs = as.character(private$run_df$used_runs)
                                  ))
          
          private$UsedEmissions <- array(dim = c(nTIMES,self$solveStates$nStates,  nRUNs),
                                         dimnames = list(
                                           time = seq(0, tmax, length.out = nTIMES),
                                           self$solveStates$asDataFrame$Abbr,
                                           RUNs = as.character(private$run_df$used_runs)
                                         ))
          
          if (needdebug) {
            debugonce(private$Function)
          }
          
          if(nRUNs > 1){
            #loop over scenarios / lhs RUNs, if any
            for (i in 1:nRUNs){
              
              # If there is one set of emissions: 
              if(!"RUN" %in% colnames(emissions) && emisFuns == 0){
                emis <- self$emissions()
                # If there are nRUNs sets of emissions or distribution functions for emissions:
              } else if(length(unique(emissions$RUN)) == nRUNs || (emisFuns != 0)){
                emis <- self$emissions(i)
                # Emissions are not of length 1 or nRUNs
              } else {
                stop("There should be 1 or nRUNs sets of emissions")
              }
              
              # possibly update dirty to create a new SB.k for uncertainty variables
              if (nVars > 0) {
                lhsruns <- private$LHSruns
                private$input_variables$Waarde <- private$LHSruns[,i]
                private$MyCore$mutateVars(private$input_variables)
                
                inputvars <- private$input_variables
                inputvars$RUN <- i
                
                vns <- unique(private$input_variables$varName)
                #update core and solve
                private$MyCore$UpdateDirty(unique(private$input_variables$varName))
                self$PrepKaasM()
                
              }
              solvedFormat <- do.call(private$Function, args = c(list(k = self$SB.k, 
                                                                      m = emis), 
                                                                 parms = list(MoreParams)))
              private$Masses[,,i] <- solvedFormat[[1]]
              private$UsedEmissions[,,i] <- solvedFormat[[2]]
              
              if(is.null(private$AllVars)){
                private$AllVars <- inputvars
              } else { 
                private$AllVars <- bind_rows(private$AllVars, inputvars) # This could maybe be done in a faster way? 
              }
              
            }
          } else { # Solve once
            emis <- self$emissions()
            solvedFormat <- do.call(private$Function, args = c(list(k = self$SB.k, 
                                                                    m = emis), 
                                                               parms = list(MoreParams)))
            dimsolved <- dim(solvedFormat[[1]])
            dimempty <- dim(private$Masses[,,1])
            private$Masses[,,1] <- solvedFormat[[1]]
            private$UsedEmissions[,,1] <- solvedFormat[[2]]
          }
        }
      },
      
      #' @description return dataframe without time and RUNs column if they 
      #' only have one unique entry
      RemoveUnusedCols = function(df){
        # Check and remove 'time' if it has only one unique value
        if (length(unique(df$time)) == 1) {
          df$time <- NULL
        }
        
        # Check and remove 'RUNs' if it has only one unique value
        if (length(unique(df$RUNs)) == 1) {
          df$RUNs <- NULL
        }
        return(df)
      },
      
      #` Function that returns the solution
      GetMasses = function() {

                # Prep and return the solution
        if (is.null(private$Masses)) {
          stop("first solve, then ask again")
        }
        SolDF <- array2DF(private$Masses)
        names(SolDF)[names(SolDF) == "Var2"] <- "Abbr"
        names(SolDF)[names(SolDF) == "Value"] <- "Mass_kg"
                
        SolDF$RUNs <- private$run_df$solver_runs[match(SolDF$RUNs, private$run_df$used_runs)]
        
        SolDF <- self$RemoveUnusedCols(SolDF)

        return(SolDF)
      },
      
      GetEmissions = function() {

        if (is.null(private$UsedEmissions)) {
          stop("first solve, then ask again")
        }
        
        EmisDF <- array2DF(private$UsedEmissions)
        names(EmisDF)[names(EmisDF) == "Var2"] <- "Abbr"
        names(EmisDF)[names(EmisDF) == "Value"] <- "Emission_kg_s"
        
        EmisDF$RUNs <- private$run_df$solver_runs[match(EmisDF$RUNs, private$run_df$used_runs)]
        
        EmisDF <- self$RemoveUnusedCols(EmisDF)
        return(EmisDF)
      },
      
      #' @description Function that returns the values for the LHS samples scaled to the distributions given by the user
      GetVarValues = function(){
        
        vars <- private$AllVars
        
        vars$RUN <- private$run_df$solver_runs[match(vars$RUN, private$run_df$used_runs)]
        
        names(vars)[names(vars) == "RUN"] <- "RUNs"
        
        return(vars)
      },
      
      #' @description Function that returns the concentration calculated from masses 
      GetConcentrations = function() {
  
        #prep and call Mass2ConcFun (Volume, Matrix, all.rhoMatrix, Fracs, Fracw)
        if (is.null(private$Masses)) {
          stop("first solve, then ask again")
        }
        if (!is.null(private$input_variables)) {
          # make sure params of private$Mass2ConcFun are not affected by dirty params
          ConcParams <- private$Mass2ConcFun$needVars
          DependVar <- private$MyCore$DependOn(ConcParams)
          if (any(unique(private$input_variables$varName) %in% c(ConcParams, DependVar))){
            stop ("concentration calculation depends on at least one of the uncertain parameters, this not implemented yet")
          }
        }

        divide <- private$Mass2ConcFun$execute()
        divide <- dplyr::left_join(private$SolveStates$asDataFrame, divide)
        solution_df <- array2DF(private$Masses) 
        names(solution_df)[names(solution_df) == "Var2"] <- "Abbr"
        solution_df <- dplyr::left_join(solution_df, divide, by="Abbr")

        solution_df$Concentration_kg_m3 <- solution_df$Value * solution_df$NewData 
        concentration_df <- solution_df[, c('time', 'RUNs', 'Abbr', 'Concentration_kg_m3')]
        
        concentration_df <- self$ConcentrationToGrams(concentration_df)
        
        concentration_df$RUNs <- private$run_df$solver_runs[match(concentration_df$RUNs, private$run_df$used_runs)]
        
        concentration_df <- self$RemoveUnusedCols(concentration_df)
        
        return(concentration_df)
      },
      
      #' @description Function that creates the appropriate concentration plot
      GetConcentrationPlot = function(scale = NULL, subcompart = NULL){
        concentration <- self$GetConcentrations() 
        
        if(is.null(scale)){
          scale <- "Regional"
          print("No scale was given to function, Regional scale is selected")
        }
        
          # Steady state deterministic
        if (identical(colnames(concentration), c("Abbr", "Concentration", "Unit"))) {
          concplot <- DetSSConcPlot(scale = scale, subcompart = subcompart)
          # Dynamic deterministic  
        } else if (identical(colnames(concentration), c("Abbr", "time",  "Concentration", "Unit"))) {
          concplot <- DetDynConcPlot(scale = scale, subcompart = subcompart)
          # Steady state probabilistic    
        } else if (identical(colnames(concentration), c("Abbr", "RUNs", "Concentration", "Unit"))) {
          concplot <- ProbSSConcPlot(scale = scale) 
          # Dynamic probabilistic  
        } else if (identical(colnames(concentration), c("Abbr", "time", "RUNs", "Concentration", "Unit"))) {
          if(is.null(subcompart)){
            subcompart <- "agriculturalsoil"
            print("No subcompart was given to function, agriculturalsoil is selected")
          }
          concplot <- ProbDynConcPlot(scale = scale, subcompart = subcompart)
        }
        
        return(concplot)
      }, 
      
      #' @description Function that creates the appropriate solution plot
      GetMassesPlot = function(scale = NULL, subcompart = NULL){
        solution <- self$GetMasses()
        
        if(is.null(scale)){
          scale <- "Regional"
          print("No scale was given to function, Regional scale is selected")
        }

          # Steady state deterministic
        if (identical(colnames(solution), c("Abbr", "Mass_kg"))) {
          solplot <- DetSSPlot(scale = scale, subcompart = subcompart)
          # Dynamic deterministic  
        } else if (identical(colnames(solution), c("time", "Abbr", "Mass_kg"))) {
          solplot <- DetDynSolPlot(scale = scale, subcompart = subcompart)
          # Steady state probabilistic    
        } else if (identical(colnames(solution), c("Abbr", "RUNs", "Mass_kg"))) {
          solplot <- ProbSSSolPlot(scale = scale) 
          # Dynamic probabilistic  
        } else if (identical(colnames(solution), c("time", "Abbr", "RUNs", "Mass_kg"))) {
          if(is.null(subcompart)){
            subcompart <- "agriculturalsoil"
            print("No subcompart was given to function, agriculturalsoil is selected")
          }
          solplot <- ProbDynSolPlot(scale = scale, subcompart = subcompart)
        }
        
        return(solplot)
      },
      
      #' @description Function that creates a mass distribution tree map for
      #' steady state solutions
      GetMassDist = function(scale = NULL){
        solution <- self$GetMasses()
        
        if(is.null(scale)){
          scale <- "Regional"
          print("No scale was given to function, Regional scale is selected")
        }
        
        # Steady state deterministic
        if (identical(colnames(solution), c("Abbr", "Mass_kg"))) {
          massdistplot <- DetSSMassDist(scale = scale)
          # Dynamic deterministic  
        } else if (identical(colnames(solution), c("time", "Abbr", "Mass_kg"))) {
          stop("No mass distribution plot available for dynamic masses")
          # Steady state probabilistic    
        } else if (identical(colnames(solution), c("Abbr", "RUNs", "Mass_kg"))) {
          massdistplot <- ProbSSMassDist(scale = scale) 
          # Dynamic probabilistic  
        } else if (identical(colnames(solution), c("time", "Abbr", "RUNs", "Mass_kg"))) {
          stop("No mass distribution plot available for dynamic masses")
        }
        
        return(massdistplot)
      },
      
      #' @description prepare kaas for matrix calculations
      #' 1 Convert relational i,j,kaas into a matrix (matrify, pivot..)
      #' including aggregation of kaas with identical (i,j)
      #' 2 The kaas are only describing removal, the where-to needs to be
      #' added to the from-masses by putting it to into the diagonal
      #' NB emission depend on order of states; if available: resort!
      #' @param kaas k's
      #' @return matrix with kaas; ready to go
      PrepKaasM = function (kaas = NULL) {

        if (is.null(kaas)) {
          #copy latest from core
          kaas <- self$myCore$kaas
        }
        if (any(kaas$k == 0.0)) {
          #or a very small value?? for solver stability?
          message(paste(table(kaas$k == 0.0)["TRUE"]), " rate constants (k values) equal to 0; removed for solver")
          kaas <- kaas[kaas$k > 0, ]
        }
        # copy, clean states (remove those without any k)
        kaas$fromIndex <- self$myCore$FindStatefrom3D(
          data.frame(
            Scale = kaas$fromScale,
            SubCompart = kaas$fromSubCompart,
            Species = kaas$fromSpecies
          )
        )
        kaas$toIndex <- self$myCore$FindStatefrom3D(
          data.frame(
            Scale = kaas$toScale,
            SubCompart = kaas$toSubCompart,
            Species = kaas$toSpecies
          )
        )
        stateInd <- sort(unique(c(kaas$fromIndex, kaas$toIndex)))
        newStates <- self$myCore$states$asDataFrame[stateInd, ]
        if (nrow(newStates) != self$myCore$states$nStates) {
          if (exists("verbose") && verbose) {
            removedStates <- do.call(paste, as.list(self$myCore$states$asDataFrame$Abbr[!self$myCore$states$asDataFrame$Abbr %in% newStates$Abbr]))
            warning(
              paste(
                self$myCore$states$nStates - nrow(newStates),
                "states without kaas, not in solver:",
                removedStates
              )
            )
          }
          private$SolveStates <- SBstates$new(newStates)
          private$SolveStates$myCore <- private$MyCore
          #redo the indices
          # copy, clean states (remove those without any k)
          kaas$fromIndex <- sapply(1:nrow(kaas), function(i) {
            which(
              newStates$Scale == kaas$fromScale[i] &
                newStates$SubCompart == kaas$fromSubCompart[i] &
                newStates$Species == kaas$fromSpecies[i]
            )
          })
          kaas$toIndex <- sapply(1:nrow(kaas), function(i) {
            which(
              newStates$Scale == kaas$toScale[i] &
                newStates$SubCompart == kaas$toSubCompart[i] &
                newStates$Species == kaas$toSpecies[i]
            )
          })
        }
        
        nrowStates <- private$SolveStates$nStates
        k2times <- as.integer(nrowStates * nrowStates)
        SB.K <- matrix(rep.int(0.0, k2times), nrow = nrowStates)
        sumkaas <-
          aggregate(k ~ fromIndex + toIndex, data = kaas, FUN = sum)
        
        for (SBi in (1:nrow(sumkaas))) {
          SB.K[sumkaas$toIndex[SBi],
               sumkaas$fromIndex[SBi]] <- sumkaas$k[SBi]
        }
        #Add the from quantities(i) to the to-states by
        #substracting the (negative) factors(i) to the diagonal
        # store the diag (== degradation and other removal processes)
        degrdiag <- diag(SB.K)
        diag(SB.K) <- 0.0 #yes, irt colSums!
        diag(SB.K) <- -degrdiag - colSums(SB.K)
        rownames(SB.K) <- newStates$Abbr
        colnames(SB.K) <- newStates$Abbr
        private$SB.K <- SB.K
        invisible(SB.K)
      },
      
      #' @description sync emissions as relational table with states into vector
      #' @param emissions named vector / 
      #' @return emissions; ready to solve for the appropriate solver
      PrepemisV = function(emis) {

        private$emissionModule <-
          EmissionModule$new(emis, private$SolveStates$asDataFrame$Abbr)
      },
      
      #' @description Function that returns the emissions for a specific RUN
      #' @param scenario run number 
      #' @return emissions for the given run
      emissions = function(scenario = NULL){
        if (is.null(private$emissionModule)) {
          stop("set emission data first, using PrepemisV()")
        }
        private$emissionModule$emissions(scenario)
      },

      PrepLHS = function(var_box_df = data.frame(), var_invFun = list(), emis_invFun = list(), nRUNs = 100){
        #checks
        # states should also be in # should be in private$SolveStates?
        solveStateAbbr <- private$SolveStates$asDataFrame
        if (!all(sapply(emis_invFun, is.function))){
          stop ("emis_invFun should be a list of functions with a single parameter")
        }
        if (!all(names(emis_invFun) %in% solveStateAbbr$Abbr)){
            stop("not all names of the emis_invFun are in states (that have k''s)")
        } 
        if (!all(sapply(var_invFun, is.function))){
          stop ("var_invFun should be a list of functions with a single parameter")
        }
        if (length(var_box_df) > 0){
          # var_box_df should contain varName and have the same length as var_invFun
          is.df.with(var_box_df, "SolverModule$PrepLHS", c("varName"))
          neededDims <- private$MyCore$fetchDims(unique(var_box_df))
          if (!all(neededDims %in% names(var_box_df))) {
            #expand from Abbr?
            var_box_df <- dplyr::left_join(var_box_df, solveStateAbbr)
          }
        }
        #only now:
        stopifnot(length(var_invFun) == nrow(var_box_df))
        
        private$input_variables <- var_box_df
        return (lhs::optimumLHS(n = nRUNs, k = length(var_invFun) + length(emis_invFun)))

      },
      
      ScaleLHS = function(lhsRUNs, var_invfun) {

        # Check if lhsRUNs is a vector and convert it to a matrix with one column if necessary
        if (is.vector(lhsRUNs)) {
          lhsRUNs <- matrix(lhsRUNs, ncol = 1)
        }
        
        # Determine the number of columns and functions
        num_columns <- ncol(lhsRUNs)
        num_functions <- length(var_invfun)
        
        # Check if there is only one function and one column
        if (num_columns == 1 && num_functions == 1) {
          # Apply the single inverse function to each element of the single column
          transformed_samples <- sapply(lhsRUNs[, 1], var_invfun[[1]])
        } else if (num_columns == num_functions) {
          # Apply each function to the corresponding column
          transformed_samples <- mapply(function(column, inv_fun) {
            sapply(column, inv_fun)
          }, as.data.frame(lhsRUNs), var_invfun, SIMPLIFY = FALSE)
          
          # Convert the list back to a matrix if needed
          transformed_samples <- do.call(cbind, transformed_samples)
        } else {
          stop("The number of columns in lhsRUNs must match the number of functions in var_invfun.")
        }
        
        return(transformed_samples)
      },
      
      #' @description diff between kaas in this and k's in OtherSB.K
      #' @param OtherSB.K the 'other' kaas
      #' @param tiny (epsilon) permitted rounding error (we might be dealing with excel/csv files ! :( )
      DiffSB.K = function(OtherSB.K, tiny = 1e-20) {
        SB.K <- self$PrepKaasM()
        #match on row/colnames?!
        rowMatch <- private$SolveStates$findState(rownames(OtherSB.K))
        colMatch <- private$SolveStates$findState(colnames(OtherSB.K))
        if (anyNA(c(rowMatch, colMatch))) {
          stop("unmatched row/col name(s) in OtherSB.K")
        }
        rowFind <-
          private$SolveStates$asDataFrame$Abbr %in% rownames(OtherSB.K)
        if (any(!rowFind)) {
          NotFound <-
            do.call(paste, as.list(private$SolveStates$asDataFrame$Abbr[!rowFind]))
          stop(paste("States missing in OtherSB.K:", NotFound))
        }
        Diff <- SB.K[rowMatch, colMatch] - OtherSB.K
        ShowDiff <- which(abs(Diff) > tiny, arr.ind = T)
        
        data.frame(
          from = private$SolveStates$asDataFrame$Abbr[ShowDiff[, 1]],
          to = private$SolveStates$asDataFrame$Abbr[ShowDiff[, 2]],
          diff = Diff[ShowDiff]
        )
      },
      
      #' @description return dataframe with concentrations
      #' in g/l, g/m3 or g/kg dw instead of kg/m3
      ConcentrationToGrams = function(Concentration_df) {
        # Fetch and filter the rhoMatrix data
        rho_data <- private$MyCore$fetchData("rhoMatrix")
        rho <- rho_data[rho_data$SubCompart == "agriculturalsoil", "rhoMatrix"]
        
        # Merge Concentration_df with MyCore states data
        Concentration_df <- merge(Concentration_df, private$MyCore$states$asDataFrame, by = "Abbr")
        
        # Create a multiplier column based on SubCompart
        Concentration_df$Multiplier <- ifelse(Concentration_df$SubCompart %in% c("air", "cloudwater"), 1000,
                                      ifelse(Concentration_df$SubCompart %in% c("river", "lake", "sea", "deepocean"), 1,
                                             ifelse(Concentration_df$SubCompart %in% c("naturalsoil", "agriculturalsoil", "othersoil", 
                                                                               "freshwatersediment", "marinesediment", "lakesediment"), 
                                                    rho * 1000, 1)))
        
        # Calculate the Concentration using the Multiplier
        Concentration_df$Concentration <- Concentration_df$Concentration_kg_m3 * Concentration_df$Multiplier
        
        # Create a Unit column based on SubCompart
        Concentration_df$Unit <- ifelse(Concentration_df$SubCompart %in% c("air", "cloudwater"), "g/m3",
                                ifelse(Concentration_df$SubCompart %in% c("river", "lake", "sea", "deepocean"), "g/L",
                                       ifelse(Concentration_df$SubCompart %in% c("naturalsoil", "agriculturalsoil", "othersoil", 
                                                                         "freshwatersediment", "marinesediment", "lakesediment"), 
                                              "g/kg dw", "kg/m3")))
        
        # Add up the concentration in air and cloudwater, name the compartment air + cloudwater
        
        # Select the desired columns
        Concentration_df <- Concentration_df[, c("Abbr", "time", "RUNs", "Concentration", "Unit")]
        return(Concentration_df)
      },
      
      #' @description return dataframe with
      #' state in three columns,
      #' time input in one or t[est]vars in separate columns,
      #' and the Mass in the Mass column
      MassesAsRelational = function(fullStates = FALSE) {
        if (is.null(private$Masses)) {
          warning("no calculation available")
          return(NULL)
        } else {
          if (fullStates){
            array2DF(private$Masses)
          } else { #append states
            arrayAsDF <- array2DF(private$Masses)
            dplyr::left_join(arrayAsDF, private$solveStates$asDataFrame)
          }
        }
      }
    ),
    active = list(
      
      #' @field needVars getter for parameters, derived from the defining function
      needVars = function(value) { #overrule
        formalArgs(private$Function)
      },
      
      #' @field states injected from States
      solveStates = function(value) { # these might differ from the core states, they are cleaned
        if (missing(value)) {
            private$SolveStates
        } else {
          stop("`$states` are set by new()", call. = FALSE)
        }
      },
      
      #' @field SB.k r.o. matrix of k's, after preparing, mostly abuse of identity for removal processes
      SB.k = function(value) {
        if (missing(value)) {
          private$SB.K
        } else {
          stop("`$SB.k` is set by PrepKaasM", call. = FALSE)
        }
      },
      
      #' @field RUNs LHS samples; TODO more generally named scenarios?
      RUNs = function(value){
        if (missing(value)){
          private$LHSruns
        } else stop("not yet possible to set RUNs, use PrepUncertain()") #or accept scenarios?
      },
      
      Input_Variables = function(value){
        private$input_variables
      }
      
    )
    
  )
