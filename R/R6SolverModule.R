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
    public = list( #
      #' @description prepare kaas for matrix calculations
      #' @param needdebug if T, the solver defining function execute in debugmode 
      execute = function(needdebug = F, emissions=NULL, solvername=NULL, ...){ 
        #browser()
        
        MoreParams <- list(...)
        
        private$SolverName <- solvername

        if (needdebug){
          debugonce(private$Function)
        }
        
        private$Solution <- do.call(private$Function, args = c(
          list(ParentModule = self),
          list(...)))
        
        solvernames <- c("EventSolver", "DynApproxSolve", "SBsolve", "UncertainSolver", "vUncertain", "UncertainDynamicSolver")
        if (private$SolverName %in% solvernames) {
          Sol <- private$Solution
          Sol
        } else {
        # Solvers (should) return a vector [states] or
        # a Matrix[states, time|run] with the mass in the state in equilibrium in the last column/row
          if (is.null(dim(private$Solution))) { #probably a vector; why is dim not 1?
            EqMass <- cbind(self$solveStates$asDataFrame, private$Solution)
          } else {
            #it's a matrix with as many rows or columns as states?
            if (! length(dim(private$Solution)) == 2 && (nrow(private$SolveStates$asDataFrame) %in% dim(private$Solution))) {
              warning("solver did not return as many rows nor cols as there are states")
              return(NULL)
            } #pick the last entry as steady state solution
            if (nrow(self$solveStates$asDataFrame) == nrow(private$Solution)) {
              private$MatrixSolutionInRows <- F
              EqMass <- cbind(self$solveStates$asDataFrame, t(as.matrix(private$Solution))[,ncol(private$Solution)])
          } else { #same amount of colums => pick the last row
            private$MatrixSolutionInRows <- T
            EqMass <- cbind(self$solveStates$asDataFrame, as.matrix(private$Solution)[nrow(private$Solution),])
          }
        }
        names(EqMass)[length(EqMass)] <- "EqMass" #last column
        EqMass 
        }
      },
      
      GetSolution = function(solvername) {
        #browser()
        if (solvername == "SB1Solve" | solvername == "SBsteady") {
          df <- as.data.frame(private$Solution)
          df <- cbind(Abbr = rownames(df), df)  # Add row names as a new column
          colnames(df) <- c("Abbr", "EqMass") 
        } else if(solvername == "UncertainSolver") {
          df <- private$Solution
        } else if(solvername == "UncertainDynamicSolver") {
          df <- private$Solution
        } else if(solvername == "DynApproxSolve") {
          df <- as.data.frame(private$Solution)
        }
        return(df)
        
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
        #browser()
        # if (!is.null(private$Emissions)){
        #   SavedEmissions <- private$Emissions
        #   private$Emissions <- NULL
        # } else { #know they should update
        #   SavedEmissions <- NULL
        # }
        
        if (is.null(kaas)) { #copy latest from core
          kaas <- self$myCore$kaas}
        if (any(kaas$k == 0.0)){ #or a very small value?? for solver stability?
          warning(paste(table(kaas$k == 0.0)["TRUE"]), " k values equal to 0; removed for solver")
          kaas <- kaas[kaas$k > 0,]
        }
        # copy, clean states (remove those without any k)
        kaas$fromIndex <- self$myCore$FindStatefrom3D(data.frame(
          Scale = kaas$fromScale,
          SubCompart = kaas$fromSubCompart,
          Species = kaas$fromSpecies
        ))
        kaas$toIndex <- self$myCore$FindStatefrom3D(data.frame(
          Scale = kaas$toScale,
          SubCompart = kaas$toSubCompart,
          Species = kaas$toSpecies
        ))
        stateInd <- sort(unique(c(kaas$fromIndex, kaas$toIndex)))
        #TODO Will the matrix be sane?
        #private$States <- SBstates$new(self$myCore$states$asDataFrame[self$myCore$states$matchi(stateInd),])
        newStates <- self$myCore$states$asDataFrame[stateInd,]
        if (nrow(newStates) != self$myCore$states$nStates) {
          if (exists("verbose") && verbose) {
            removedStates <- do.call(paste, as.list(
                        self$myCore$states$asDataFrame$Abbr[!self$myCore$states$asDataFrame$Abbr %in% newStates$Abbr]
            ))
            warning(paste(self$myCore$states$nStates - nrow(newStates),"states without kaas, not in solver:", removedStates))
            #remove states without state in columns OR rows => matrix is singular 
          }
          private$SolveStates <- SBstates$new(newStates) 
          private$SolveStates$myCore <- private$MyCore
          #redo the indices
          # copy, clean states (remove those without any k)
          kaas$fromIndex <- sapply(1:nrow(kaas), function(i){
            which(newStates$Scale == kaas$fromScale[i] &
                  newStates$SubCompart == kaas$fromSubCompart[i] &
                  newStates$Species == kaas$fromSpecies[i])
          })
          kaas$toIndex <- sapply(1:nrow(kaas), function(i){
            which(newStates$Scale == kaas$toScale[i] &
                    newStates$SubCompart == kaas$toSubCompart[i] &
                    newStates$Species == kaas$toSpecies[i])
          })
        }  
        
        #Reinstate(s) emissions
        # if (!is.null(SavedEmissions)){
        #   OldEmissions <- data.frame(
        #     Abbr = names(SavedEmissions)[!is.na(names(SavedEmissions))],
        #     Emis = SavedEmissions[!is.na(names(SavedEmissions))]
        #   ) 
        #   self$PrepemisV(OldEmissions)
        # }
        
        nrowStates <- self$solveStates$nStates
        k2times <- as.integer(nrowStates*nrowStates)
        SB.K <- matrix(rep.int(0.0, k2times), nrow = nrowStates)
        sumkaas <- aggregate(k ~ fromIndex + toIndex, data = kaas, FUN = sum)

        for(SBi in (1:nrow(sumkaas))){
          SB.K[sumkaas$toIndex[SBi],
               sumkaas$fromIndex[SBi]] <- sumkaas$k[SBi] 
        }
        #Add the from quantities(i) to the to-states by
        #substracting the (negative) factors(i) to the diagonal
        # store the diag (== degradation and other removal processes)
        degrdiag <- diag(SB.K)
        diag(SB.K) <- 0.0 #yes, irt colSums!
        diag(SB.K) <- - degrdiag - colSums(SB.K)
        rownames(SB.K) <- newStates$Abbr
        colnames(SB.K) <- newStates$Abbr
        private$SB.K <- SB.K
        invisible(SB.K)
      },

      #' @description sync emissions as relational table with states into vector 
      #' @param emissions 
      #' @return vector of functions(t) with emissions; ready to solve
      PrepemisV = function (emissions = NULL, solvername = private$SolverName, ...) { #for backward compatibly
        #browser()

        private$Emission <- EmissionModule$new(self, emissions, solvername, private$SB.K, ...)

        emis <- private$Emission$CleanEmissions()
        private$Emissions <- private$Emission$CleanEmissions()

      },
      
      PrepUncertain = function(input) {
        #browser()
                # Colnames that should be in the df
        cn <- c("varName", "Scale", "SubCompart", "data")
        
        if(!all(cn %in% names(input))) {
          stop("Column name(s) incorrect. The tibble should contain columns with the following names: 'varname', 'scale', 'subcompart' and 'data'.")
        }
        
        row_counts <- input %>%
          pull(data) %>%
          map_int(nrow)
        
        unique_rc <- length(unique(row_counts)) 
        
        if(unique_rc != 1) {
          stop("Not all variables have the same number of samples.")
        }
        
        private$lUncertainInput <- input
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
      
      EmisBoxODE = function(t, m, parms){
        e_t <- parms$e(t)
        dm <- with(parms, K %*% m + e_t)
        return(list(dm))
      },
      
      EventODE = function(t, m, parms){
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
        rowMatch <- self$solveStates$findState(rownames(OtherSB.K))
        colMatch <- self$solveStates$findState(colnames(OtherSB.K))
        if (anyNA(c(rowMatch, colMatch))) {
          #browser()
          #unique(is.NA(rowMeans))
          stop("unmatched row/col name(s) in OtherSB.K")
        }
        rowFind <- self$solveStates$asDataFrame$Abbr %in% rownames(OtherSB.K)
        if (any(!rowFind)) {
          NotFound <- do.call(paste, as.list(self$solveStates$asDataFrame$Abbr[!rowFind]))
          stop(paste("States missing in OtherSB.K:", NotFound))
        }
        Diff <- SB.K[rowMatch, colMatch] - OtherSB.K
        ShowDiff <- which(abs(Diff) > tiny, arr.ind = T)
        
        data.frame(
          from = self$solveStates$asDataFrame$Abbr[ShowDiff[,1]],
          to = self$solveStates$asDataFrame$Abbr[ShowDiff[,2]],
          diff = Diff[ShowDiff]
        )
      },
      #' @description return dataframe with 
      #' state in three columns, 
      #' time input in one or t[est]vars in separate columns, 
      #' and the Mass in the Mass column
      SolutionAsRelational = function(...){
        if (is.null(self$SBtime.tvars) && is.null(self$vnamesDistSD)) {
          warning("no calculation available")
          return(NULL)
        }
        if (!is.null(self$SBtime.tvars)) {
          if (private$MatrixSolutionInRows) {
            ret2Blong <- cbind(private$SolveStates$asDataFrame[,The3D], 
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
            if(length(Params) == 0) { #default uncertainty; return mean and sd.
              #now we know the structure / names
              #exclude base calculation from the stats
              exclNr <- which(rownames(private$Solution) == "base")
              clMeans <- colMeans(private$Solution[-exclNr,])
              clSD <- apply(private$Solution[-exclNr,], 2, sd)
              pval <- list()
              for (px in c(0.05, 0.25, 0.75, 0.95)) {
                percname <- paste("p", round(100*px), sep = "")
                pval[[percname]] <- apply(private$Solution[-exclNr,], 2, quantile, probs = px)
              }
              ret <- do.call(rbind, pval)
              ret <- as.data.frame(t(rbind(clMeans, clSD, ret)))
              return(cbind(self$solveStates$asDataFrame[,The3D], ret))
              
            } else {
              #is there a formula "terms"
              if ("terms" %in% names(Params)) {
                FUN <- Params$FUN
                if (is.null(FUN)) {
                  warning("no FUN found, no aggregation")
                  #browser()
                  leftRightHS <- private$interprFormula(Params$terms)
                  #Add Name to attribute name, if needed
                  fNasName <- sapply(fAsList, function(x) {
                    #browser()
                    ifelse (any(x %in% The3D), paste(x, "Name", sep = ""), x)
                  })
                  
                  ParsAndDims <- sapply(fAsList, function(x)
                    self$SolveStates$findDim(x))
                  
                } else {
                  #call and return aggregate
                  aggregate(x = Params$terms, data = private$Solution, FUN = FUN)
                }
              }
              if (length(Params) == 1 && is.character(Params[[1]])) {
                VarName <- Params[[1]]
                if (! VarName %in% self$vnamesDistSD$vnames) {
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
          if (is.null(private$SolveStates)) {
            return(self$myCore$states)
          } else {
            private$SolveStates
          }
        } else {
          stop("`$states` are set by new()", call. = FALSE)
        }
      },
      
      #' @field emissions vector of emissions
      emissions = function(value) { 
        if (missing(value)) {
          private$Emissions
        } else {
          
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
      
      UncertainInput = function(value){
        private$lUncertainInput
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
    ),
    
    private = list(
      
      #' @description return list with vectors of 
      #' subcompartments, variables and filters from subcompartments ~ variables | filters
      interprFormula = function(TheCall){ 
        
        plusflat <- function(lisTree){
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
      ,
      NeedVars = function(){ #overrule
        NULL
      },
      Solution = NULL,
      ST = NULL,
      SolveStates = NULL,
      Emissions = NULL,
      emis = NULL, 
      SB.K = NULL,
      MatrixSolutionInRows = NULL,
      lSBtime.tvars = NULL,
      lvnamesDistSD = NULL,
      Emission = NULL,
      SolverName = NULL,
      lUncertainInput = NULL
    )
  )
