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
      execute = function(needdebug = F){ 
        if (needdebug){
          debugonce(private$Function)
        }
        private$Solution <- do.call(private$Function, args = c(
          list(ParentModule = self),
          private$MoreParams))
        # Solvers (should) return a vector [states] or
        # a Matrix[states, time]
        lastCol <- ncol(as.matrix(private$Solution))
        EqMass <- cbind(private$States$asDataFrame, as.matrix(private$Solution)[,lastCol])
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

        # copy, clean states (remove those without any k)
        stateInd <- unique(c(kaas$i, kaas$j))
        private$States <- SBstates$new(self$myCore$states$asDataFrame[self$myCore$states$matchi(stateInd),])
        if (private$States$nStates != self$myCore$nStates && exists("verbose") && verbose) {
          warning(paste(self$myCore$nStates - private$States$nStates,"without kaas, not in solver"))
        }
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
        # store the diag (== degradation)
        degrdiag <- diag(SB.K)
        diag(SB.K) <- 0.0 #yes, irt colSums!
        diag(SB.K) <- - degrdiag - colSums(SB.K)
        private$SB.K <- SB.K
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
      }
    ),
    
    private = list(
      NeedVars = function(){ #overrule
        NULL
      },
      Solution = NULL,
      States = NULL,
      Emissions = NULL,
      SB.K = NULL
    )
  )
