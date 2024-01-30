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
      execute = function(needdebug = F, ...){ 
        if (needdebug){
          debugonce(private$Function)
        }
        private$Solution <- do.call(private$Function, args = c(
          list(ParentModule = self),
          list(...)))
        # old ?
        # private$Solution <- do.call(private$Function, args = c(
        #   list(ParentModule = self),
        #   private$MoreParams))
        # Solvers (should) return a vector [states] or
        # a Matrix[states, time|run] with the mass in the state in equilibrium in the last column/row
        if (is.null(dim(private$Solution))) { #probably a vector; why is dim not 1?
          EqMass <- cbind(private$States$asDataFrame, private$Solution)
        } else {
          #it's a matrix with as many rows or columns as states?
          if (! length(dim(private$Solution)) == 2 && (nrow(private$States$asDataFrame) %in% dim(private$Solution))) {
            warning("solver did not return as many rows nor cols as there are states")
            return(NULL)
          } #pick the last entry as steady state solution
          if (nrow(private$States$asDataFrame) == nrow(private$Solution)) {
            private$MatrixSolutionInRows <- F
            EqMass <- cbind(private$States$asDataFrame, t(as.matrix(private$Solution))[,ncol(private$Solution)])
          } else { #same amount of colums => pick the last row
            private$MatrixSolutionInRows <- T
            EqMass <- cbind(private$States$asDataFrame, as.matrix(private$Solution)[nrow(private$Solution),])
          }
        }
        names(EqMass)[length(EqMass)] <- "EqMass" #last column
        EqMass 
      },
      
      #' @description prepare kaas for matrix calculations
      #' 1 Convert relational i,j,kaas into a matrix (matrify, pivot..)
      #' including addition of kaas with identical (i,j)
      #' 2 The kaas are only describing removal, the where-to needs to be
      #' added to the from-masses by putting it to into the diagonal
      #' @param kaas k's
      #' @return matrix with kaas; ready to go
      PrepKaasM = function (kaas = NULL) {
        if (is.null(kaas)) { #copy latest from core
          kaas <- self$myCore$kaas}

        if (any(kaas$k == 0.0)){ #or a very small value?? for solver stability?
          warning(paste(table(kaas$k == 0.0)["TRUE"]), " k values equal to 0; removed for solver")
          kaas <- kaas[kaas$k > 0,]
        }
        # copy, clean states (remove those without any k)
        stateInd <- unique(c(kaas$i, kaas$j))
        #private$States <- SBstates$new(self$myCore$states$asDataFrame[self$myCore$states$matchi(stateInd),])
        newStates <- self$myCore$states$asDataFrame[self$myCore$states$matchi(stateInd),]
        if (nrow(newStates) != self$myCore$states$nStates && exists("verbose") && verbose) {
          warning(paste(self$myCore$states$nStates - nrow(newStates),"states without kaas, not in solver"))
        }
        #remove states without state in columns OR rows => matrix is singular 
        #TODO  dynamic calculation is still possible! if emission is a source??
        stateInd <- unique(kaas$i) 
        stateIndBoth <- stateInd[stateInd %in% unique(kaas$j)]
        if (length(stateInd) != length(stateIndBoth)) {
          warning(paste(length(stateInd) - length(stateIndBoth), " states with no sources OR no sinks"))
          newStates <- newStates[newStates$i %in% stateIndBoth,]
        }
        private$States <- SBstates$new(newStates)
        nrowStates <- private$States$nStates
        k2times <- as.integer(nrowStates*nrowStates)
        SB.K <- matrix(rep.int(0.0, k2times), nrow = nrowStates)
        sumkaas <- aggregate(k ~ i + j, data = kaas, FUN = sum)
        for(SBi in (1:nrow(sumkaas))){
          SB.K[private$States$matchi(sumkaas$j[SBi]),
               private$States$matchi(sumkaas$i[SBi])] <- sumkaas$k[SBi] 
        }
        #Add the from quantities(i) to the to-states by
        #substracting the (negative) factors(i) to the diagonal
        # store the diag (== degradation and other removal processes)
        degrdiag <- diag(SB.K)
        diag(SB.K) <- 0.0 #yes, irt colSums!
        diag(SB.K) <- - degrdiag - colSums(SB.K)
        private$SB.K <- SB.K
        invisible(SB.K)
      },

      #' @description sync emissions as relational table with states into vector 
      #' @param emissions 
      #' @return matrix with kaas; ready to go
      PrepemisV = function (emissions) {
        if (is.null(private$States)) {
          stop("kaas should be set first, using PrepKaasM(), to set the clean states")
        }
        if (!("data.frame" %in% class(emissions)))
          stop("emssions are expected in a data.frame-like class", call. = FALSE)
        if (!all(c("Abbr","Emis") %in% names(emissions)))
          stop("emissions should contain 'Abbr','Emis'. Note the capitals.", call. = FALSE)
        
        #private$Emissions <- value[,c('Abbr','Emis')]
        
        vEmis <- rep(0.0, length.out = private$States$nStates)
        names(vEmis)[private$States$findState(emissions$Abbr)] <- emissions$Abbr
        vEmis[private$States$findState(emissions$Abbr)] <- emissions$Emis
        # from kg/yr to Mol/s
        Molweight <- self$myCore$fetchData("MW")
        private$Emissions <- vEmis * 1000000 / Molweight / (3600*24*365) #t/an -> mol/s
        private$Emissions
      },
      
      #' @description basic ODE function for simplebox; the function for the ode-call; 
      #' see desolve/rootsolve packages
      #' @param t time (vector ?)
      #' @param m  (i) = initial mass
      #' @param parms = kaas, emissions
      #' @returns dm (i) = change in mass as list
      SimpleBoxODE = function(t, m, parms) {
        dm <- with(parms, K %*% m + e)
        return(list(dm))
      },
      
      #' @description diff between kaas in this and k's in OtherSB.K
      #' @param OtherSB.K the 'other' kaas
      #' @param tiny (epsilon) permitted rounding error (we might be dealing with excel ! :( )
      DiffSB.K = function(OtherSB.K, tiny = 1e-20) {
        SB.K <- self$PrepKaasM()
        #match on row/colnames?!
        rowMatch <- private$States$findState(rownames(OtherSB.K))
        colMatch <- private$States$findState(colnames(OtherSB.K))
        if (anyNA(c(rowMatch, colMatch))) {
          stop("unmatched row/col name(s) in OtherSB.K")
        }
        Diff <- SB.K[rowMatch, colMatch] - OtherSB.K
        ShowDiff <- which(abs(Diff) > tiny, arr.ind = T)
        
        data.frame(
          from = private$States$asDataFrame$Abbr[ShowDiff[,1]],
          to = private$States$asDataFrame$Abbr[ShowDiff[,2]],
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
            ret2Blong <- cbind(private$States$asDataFrame[,The3D], 
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
              browser()
              #now we know the structure / names
              #exclude base calculation from the stats
              exclNr <- which(rownames(private$Solution) == "base")
              clMeans <- colMeans(private$Solution[-exclNr,])
              clSD <- apply(private$Solution[-exclNr,], 2, sd)
              pval <- list()
              for (px in c(0.05, 0.25, 0.75, 0.95)) {
                percname <- paste("p", round(100*px), sep = "")
                pval[[percname]] <- apply(private$Solution[-exclNr,], 2, quantile, props = px)
              }
              ret <- do.call(rbind, pval)
              ret <- rbind(rwMeans, clSD, ret)
            } else {
              #is there a formula "terms"
              if ("terms" %in% names(Params)) {
                FUN <- Params$FUN
                if (is.null(FUN)) {
                  warning("no FUN found, no aggregation")
                  browser()
                  leftRightHS <- private$interprFormula(Params$terms)
                  #Add Name to attribute name, if needed
                  fNasName <- sapply(fAsList, function(x) {
                    browser()
                    ifelse (any(x %in% The3D), paste(x, "Name", sep = ""), x)
                  })
                  
                  ParsAndDims <- sapply(fAsList, function(x)
                    self$states$findDim(x))
                  
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
      states = function(value) { # these might differ from the core states, they are cleaned
        if (missing(value)) {
          private$States
        } else {
          stop("`$states` are set by new()", call. = FALSE)
        }
      },
      
      #' @field emissions guess what?
      emissions = function(value) { 
        if (missing(value)) {
          private$Emissions
        } else {
          stop("`$emissions` are set by PrepemisV", call. = FALSE)
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
        browser()
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
      States = NULL,
      Emissions = NULL,
      SB.K = NULL,
      MatrixSolutionInRows = NULL,
      lSBtime.tvars = NULL,
      lvnamesDistSD = NULL
    )
  )
