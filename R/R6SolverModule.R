#' @title SolverModule
#' @description
#' wrapper for your solver function and collection of general methods needed for most solvers,
#' like preparing kaas matrix and matching emissions vector and
#' the general ODE for simplebox
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
      Solution = NULL,
      SolveStates = NULL,
      EmissionModule = NULL,
      SB.K = NULL,
      MatrixSolutionInRows = NULL,
      lSBtime.tvars = NULL,
      lvnamesDistSD = NULL,
      #SolverOutputTp = NULL,
      input_variables = NULL,
      Mass2ConcFun = NULL,
      LHSruns = NULL,
      
      #' @description helper function for Make_inv_unif01
      triangular_cdf_inv = function(u, # LH scaling factor
                                    a, # Minimum
                                    b, # Maximum
                                    c) { # Peak value
        ifelse(u < (c-a)/(b-a),
               a + sqrt(u * (b-a) * (c-a)),
               b - sqrt((1-u) * (b-a) * (b-c)))
      },
      
      
      #' @description return list with vectors of
      #' subcompartments, variables and filters from subcompartments ~ variables | filters
      interprFormula = function(TheCall) {
        plusflat <- function(lisTree) {
          if (is.call(lisTree)) {
            return(c(plusflat[2], plusflat[3]))
          } else {
            return(structure(listree)) #remove formule attributes
          }
        }
        #browser()
        if (is.call(TheCall[1]) && TheCall[[1]] == "~") {
          hasTilde <- T
          lhs <- plusflat(TheCall[2])
        } else {
          hasTilde <- F
          rhs <- plusflat(TheCall[3])
        }
        list(lhs = lhs, rhs = rhs)
      }
      
      
    ),
    
    public = list(
      initialize = function(TheCore, exeFunction, ...) {
        # Call the parent class's initialize method; inits also private$MyName
        super$initialize(TheCore, exeFunction, ...)
        
        # a SB variable outside of the ModuleList and nodelist of the core
        private$Mass2ConcFun <-
          VariableModule$new(TheCore, "Mass2Conc", IsVectorised = T, AggrBy = NULL, AggrFun = NULL)
      },
      #' @description Solve matrix / emsissions
      #' @param needdebug if T, the solver defining function execute in debugmode
      execute = function(needdebug = F,
                         emissions = NULL,
                         ...) {

        MoreParams <- list(...)
        
        if (is.null(private$SB.K)) {
          stop("run PrepKaasM() first")
        }
        
        if (!is.null(emissions)) {
          self$PrepemisV(emissions)
        } else {
          if (is.null(private$EmissionModule)){
            stop ("emissions are missing")
          }
        }
        
        if (needdebug) {
          debugonce(private$Function)
        }
        
        private$Solution <- do.call(private$Function, args = c(list(ParentModule = self),
                                                               list(...)))
        
        # if ("SteadyStateMass" %in% names(private$Solution)){
        #   private$Input_Variables = private$Solution$Input_Variables
        #   private$Input_Emission = private$Solution$Input_Emission
        #   private$Solution = private$Solution$SteadyStateMass
        #   private$SolverOutputTp <- "SteadyStateMass"
        #   return(list(Input_Variables = private$Input_Variables, 
        #               Input_Emission = private$Input_Emission, 
        #               SteadyStateMass = private$Solution))
        # }
        # 
        # # check and derive SolverOutputTp from solveroutput
        # if (is.null(dim(private$Solution))) {
        #   private$SolverOutputTp <- "eqVector"
        #   private$Solution <-
        #     cbind(private$SolveStates$asDataFrame, private$Solution)
        #   return(private$Solution)
        # }
        # 
        # # it's a matrix with as many rows or columns as states?
        # if (!tibble::is_tibble(private$SolverOutputTp) & length(dim(private$Solution)) == 2) {
        #   private$SolverOutputTp <- "massMatrix"
        #   # attach states
        #   attributes(private$SolverOutputTp) <-
        #     private$SolveStates$asDataFrame
        #   return(private$Solution)
        # }
        # 
        # browser()
        # # a nested tibble ?
        # if (tibble::is_tibble(private$SolverOutputTp) &&
        #     any(sapply(x, is.list))) {
        #   private$SolverOutputTp <- "nestedTibble"
        #   # attach states
        #   attributes(private$SolverOutputTp) <-
        #     private$SolveStates$asDataFrame
        #   return(private$Solution)
        # }
        # 
        # #TODO 4e type
        
        
      },
      
      GetSolution = function(solvername) {
        return(private$Solution)
      },
      
      GetConcentrations = function() {
        #prep and call Mass2ConcFun (EqMass (as squeeze), Volume, Matrix, all.rhoMatrix, Fracs, Fracw)
        if (is.null(private$Solution)) {
          stop("first solve, then ask again")
        }
        # switch (
        #   private$SolverOutputTp,
        #   "eqVector" = {
        #     debugonce(private$Mass2ConcFun)
        #     prepped <- private$Mass2ConcFun(squeezeVar = private$Solution)
        #     browser()
        #   },
        #   "massMatrix" = {
        #     browser()
        #     debugonce(private$Mass2ConcFun)
        #     prepped <- private$Mass2ConcFun(squeezeVar = private$Solution)
        #   },
        #   "nestedTibble" = {
        #     browser()
        #     debugonce(private$Mass2ConcFun)
        #     prepped <- private$Mass2ConcFun(squeezeVar = private$Solution)
        #   }
        # )
        
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
          warning(paste(table(kaas$k == 0.0)["TRUE"]), " k values equal to 0; removed for solver")
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
        #TODO Will the matrix be sane?
        #private$States <- SBstates$new(self$myCore$states$asDataFrame[self$myCore$states$matchi(stateInd),])
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
            #remove states without state in columns OR rows => matrix is singular
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
      PrepemisV = function (emissions) {
        private$EmissionModule <-
          EmissionModule$new(emissions)
      },
      
      emissions = function(){
        if (is.null(private$EmissionModule)) {
          stop("set emission data first, using PrepemisV()")
        }
        private$EmissionModule$emissions(private$SolveStates)
      },
      
      emissionFunctions = function(){
        if (is.null(private$EmissionModule)) {
          stop("set emission data first, using PrepemisV()")
        }
        private$EmissionModule$emissionFunctions(private$SolveStates)
      },

      #create a function for transformation of lhs range (0-1) to actual variable range (inverse of the 0-1 cdf)
      Make_inv_unif01 = function(fun_type = "triangular", pars) {
        if (!fun_type %in% c("triangular", "norm", "unif")) {
          stop("! fun_type %in% c('triangular', 'norm', 'unif')")
        }
        if (fun_type == "triangular") {
          if (!(inherits(pars, "list") && length(pars) == 3)) {
            stop(
              "the triangular is created using a list of three parameters, a = minimum, b = maximum, c = peak")
          }
          a <- pars[["a"]]
          b <- pars[["b"]]
          c <- pars[["c"]]
          return(function(x) {
            private$triangular_cdf_inv(x, a, b, c)
          })
        }
        if (fun_type == "norm") {
          if (!(inherits(pars, "list")) && length(pars) == 2) {
            stop("the norm is created using a list of two parameters, mean = mean, sd = maximum, c = peak")
          }
          mu <- pars[["mean"]]
          sig <- pars[["sigma"]]
          return(function(x) {
            qnorm(p = x, mean = mu, sd = sig)
          })
        }
        if (fun_type == "unif") {
          if (!(inherits(pars, "list")) && length(pars) == 2) {
            stop("the unif is created using a list of two parameters, min = minimum, max = maximum")
          }
          minx <- pars[["min"]]
          maxx <- pars[["max"]]
          return(function(x) {
            minx + x * (maxx - minx)
          })
        }
        
      },
        
      PrepUncertain = function(varname_states, emis_states, var_invFun, emis_invFun) {

        #browser()
        # varName & Dims of all vars present
        is.df.with(varname_states, callingName = "SolverModule$PrepUncertain", 
                   mustHavecols = c('varName'))
        allVarDims <- unique(varname_states$varName)
        
        PrepemisV(emis_states)
        
        row_counts <- input %>%
          pull(data) %>%
          map_int(nrow)
        
        unique_rc <- length(unique(row_counts))
        
        if (unique_rc != 1) {
          stop("Not all variables have the same number of samples.")
        }
        
        sample_df <- input |> 
            mutate(nRUNs = map_int(data, nrow)) |> 
            mutate(
              data = map(data, ~ .x |> 
                           mutate(RUN = 1:unique(nRUNs)))
            ) |> select(-nRUNs)
        
        private$input_variables <- sample_df

      },
      
      PrepLHS = function(emis_scen = NULL, var_uncertain = NULL, nRUNs = 100){
        #checks
        if (!is.null(emis_scen)){ 
          # Either Abbr or all The3D should be in 
          if ("Abbr" %in% names(emis_scen)){
            test(is.df.with(emis_scen, "SolverModule.PrepLHS", c("timed", "emis", "Abbr")))
          } else {
            test(is.df.with(emis_scen, "SolverModule.PrepLHS", c("timed", "emis", The3D)))
            #has_scenario <- "scenario" %in% names(emis_scen)
          }
        } 
        if (!is.null(var_uncertain)){
          # Either Abbr or variable-needed The3D should be in 
          if ("Abbr" %in% names(var_uncertain)){
            test(is.df.with(var_uncertain, "SolverModule.PrepLHS", c("varName", "waarde", "Abbr")))
          } else {
            test(is.df.with(var_uncertain, "SolverModule.PrepLHS", c("varName", "waarde", The3D)))
            #has_scenario <- "scenario" %in% names(var_uncertain)
          }
        }
        # states should also be in # should be in private$SolveStates?
        private$LHSruns <- lhs::optimumLHS()

      },
      
      #' @description basic ODE function for simplebox; the function for the ode-call;
      #' see desolve/rootsolve packages
      #' @param t time (vector ?)
      #' @param m  (i) = initial mass
      #' @param parms = kaas, emissions
      #' @returns dm (i) = change in mass as list
      SimpleBoxODE = function(t, m, parms) {
        dm <- with(parms, K %*% m + e)
        return(list(dm, signal = parms$e))
      },
      
      EmisBoxODE = function(t, m, parms) {
        e_t <- parms$e(t)
        dm <- with(parms, K %*% m + e_t)
        return(list(dm))
      },
      
      EventODE = function(t, m, parms) {
        with(as.list(c(parms, m)), {
          dm <- K %*% m
          list(c(dm))
        })
      },
      
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
          #browser()
          #unique(is.NA(rowMeans))
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
      #' @description return dataframe with
      #' state in three columns,
      #' time input in one or t[est]vars in separate columns,
      #' and the Mass in the Mass column
      SolutionAsRelational = function(...) {
        if (is.null(self$SBtime.tvars) && is.null(self$vnamesDistSD)) {
          warning("no calculation available")
          return(NULL)
        }
        if (!is.null(self$SBtime.tvars)) {
          if (private$MatrixSolutionInRows) {
            ret2Blong <- cbind(private$SolveStates$asDataFrame[, The3D],
                               t(as.matrix(private$Solution)))
          }
          ret = tidyr::pivot_longer(ret2Blong,
                                    names(ret2Blong)[!names(ret2Blong) %in% The3D],
                                    names_to = "SBtime",
                                    values_to = "Mass")
          #replace the generated columnnames from SBtime into actual times
          TheTimes <- unlist(self$SBtime.tvars)
          ret$SBtime <- TheTimes[as.numeric(ret$SBtime)]
          return(ret)
        } else {
          if (!is.null(self$vnamesDistSD)) {
            Params <- list(...)
            if (length(Params) == 0) {
              #default uncertainty; return mean and sd.
              #now we know the structure / names
              #exclude base calculation from the stats
              exclNr <- which(rownames(private$Solution) == "base")
              clMeans <- colMeans(private$Solution[-exclNr, ])
              clSD <- apply(private$Solution[-exclNr, ], 2, sd)
              pval <- list()
              for (px in c(0.05, 0.25, 0.75, 0.95)) {
                percname <- paste("p", round(100 * px), sep = "")
                pval[[percname]] <-
                  apply(private$Solution[-exclNr, ], 2, quantile, probs = px)
              }
              ret <- do.call(rbind, pval)
              ret <- as.data.frame(t(rbind(clMeans, clSD, ret)))
              return(cbind(private$SolveStates$asDataFrame[, The3D], ret))
              
            } else {
              #is there a formula "terms"
              if ("terms" %in% names(Params)) {
                FUN <- Params$FUN
                if (is.null(FUN)) {
                  warning("no FUN found, no aggregation")
                  #browser()
                  leftRightHS <-
                    private$interprFormula(Params$terms)
                  #Add Name to attribute name, if needed
                  fNasName <- sapply(fAsList, function(x) {
                    #browser()
                    ifelse (any(x %in% The3D), paste(x, "Name", sep = ""), x)
                  })
                  
                  ParsAndDims <- sapply(fAsList, function(x)
                    private$SolveStates$findDim(x))
                  
                } else {
                  #call and return aggregate
                  aggregate(
                    x = Params$terms,
                    data = private$Solution,
                    FUN = FUN
                  )
                }
              }
              if (length(Params) == 1 &&
                  is.character(Params[[1]])) {
                VarName <- Params[[1]]
                if (!VarName %in% self$vnamesDistSD$vnames) {
                  warning(paste(VarName, "not in analyses"))
                  return(data.frame(NA))
                }
                #Yes, we can
                
              } else {
                warning("expected a variable name")
              }
              
            }
            
          }
        } # TODO else t[est]vars : uncertainty / sensitivity solvers etc.
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
        } else stop("not yet possible to set RUNs, use PrepUncertain()")
      },
      
      Input_Variables = function(value){
        private$input_variables
      },
      
      SBtime.tvars = function(value) {
        if (missing(value)) {
          private$lSBtime.tvars
        } else {
          private$lSBtime.tvars <- value
          private$lvnamesDistSD <- NULL
        }
      },
      
      vnamesDistSD = function(value){
        if (missing(value)) {
          private$lvnamesDistSD
        } else {
          private$lvnamesDistSD <- value
          private$lSBtime.tvars <- NULL
        }
      }
    )
    
  )
