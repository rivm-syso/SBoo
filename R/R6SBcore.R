#' @title SBcore
#' @description the centre object for running simple box with an R script.
#' @import R6
#' @export
SBcore <- R6::R6Class("SBcore",
  public = list(
    #' @description init
    #' @param NewstateModule The stateModule with its States and standard input data
    initialize = function(NewstateModule){
      private$SB4Ndata <- NewstateModule$SB4N.data
      private$States <- SBstates$new(NewstateModule$states) 
      private$substance <- NewstateModule$substance
      private$ModuleList <- list()
    },
    #' @description add a process to the calculations
    #' @param ProcessFunction The name (as character) of the process defining function 
    #' see 
    NewProcess = function(ProcessFunction){
      #existing function?
      stopifnot("function" %in% class(match.fun(ProcessFunction)))
      #delete, if it exists (being updated)
      if (ProcessFunction %in% names(private$ModuleList)){
        private$ModuleList <- private$ModuleList[private$ModuleList != ProcessFunction]
        private$nodeList <- private$nodeList[private$nodeList$Calc != ProcessFunction,]
      }
      aNewProcessModule <- ProcessModule$new(self,ProcessFunction)
      private$storeNodes(aNewProcessModule)
      private$ModuleList[[ProcessFunction]] <- aNewProcessModule
      invisible(aNewProcessModule)
    },
    
    #' @description remove a set of SB items from the calculations, for later reprocessing
    #' @param VarFunctions The name (as character) of the process defining function 2be postponed
    #' @param FlowFunctions The name (as character) of the process defining function 2be postponed
    #' @param ProcesFunctions The name (as character) of the process defining function 2be postponed
    #' @details The SBitems will be removed from the internal structures, but the names of their defining functions 
    #' will be stored locally. The private method DoPostponed will add and execute them. (This is automised, 
    #' because it will be needed in sensitivity analyses etc.)
    PostponeVarProcess = function(VarFunctions = NULL, FlowFunctions = NULL, ProcesFunctions) {
      private$l_postPoneList <- list(VarFunctions, FlowFunctions, ProcesFunctions)
    },
    
    #' @description add a flow to the calculations
    #' @param FlowFunction The name (as character) of the flow defining function 
    #' see 
    #' @param WithProcess The name (as character) of the process related to this function
    #' This is normally the k_Advection process; created if not already present
    NewFlow = function(FlowFunction, WithProcess = "k_Advection"){
      #existing function?
      stopifnot("function" %in% class(match.fun(FlowFunction)))
      
      #delete FlowModule, if it exists already (being updated)
      if (FlowFunction %in% names(private$ModuleList)){
        private$ModuleList <- private$ModuleList[private$ModuleList != FlowFunction]
        private$nodeList <- private$nodeList[private$nodeList$Calc != FlowFunction,]
      }
      aNewFlowModule <- FlowModule$new(self, FlowFunction, WithProcess)
      private$storeNodes(aNewFlowModule)

      # book-keeping for the WithProcess
      stopifnot("function" %in% class(match.fun(WithProcess)))
      # is it an advection process, i.e. proper parameter(s)
      stopifnot("flow" %in% formalArgs(WithProcess))
      # make the process, if not present
      if (WithProcess %in% names(private$ModuleList)) {
        nodeProcess <- private$ModuleList[[WithProcess]]
      } else {
        nodeProcess <- self$NewProcess(WithProcess)
      }
      if(length(nodeProcess$WithFlow)==1 && is.na(nodeProcess$WithFlow)){
        nodeProcess$WithFlow <- FlowFunction
      } else {
        nodeProcess$WithFlow <- c(nodeProcess$WithFlow, FlowFunction) 
      }
      # update needVars and private$nodeList for WithProcess; 
      # replace flow (the parameter in the function) with all individual flows, for DAG
      OrigVars <- formalArgs(WithProcess)
      nodeProcess$needVars <- c(OrigVars[OrigVars != "flow"], nodeProcess$WithFlow)
      #possibly replace existing Params "flow", or add the new flowname
      indProcess <- which(WithProcess == private$nodeList$Calc)
      matchWithFlow <- match("flow", private$nodeList$Params[indProcess])
      if (!is.na(matchWithFlow) && matchWithFlow > 0) {
        private$nodeList$Params[indProcess[matchWithFlow]] <- FlowFunction
      } else {
        private$nodeList[nrow(private$nodeList)+1,] <- c(WithProcess, FlowFunction, "Process")
      }
      private$ModuleList[[FlowFunction]] <- aNewFlowModule
      invisible(aNewFlowModule)
      
    },
    
    #' @description add a SBVariable to the calculations
    #' @param VariableFunction The name (as character) of the SBvariable defining function 
    #' see 
    #' @param AggrBy exceptional parameter to aggregate over a dim (1 of The3D)
    #' @param AggrFun FUN to use in aggregation over a dim (1 of The3D)
    NewCalcVariable = function(VariableFunction, AggrBy = NA, AggrFun = NA){
      #delete, if it exists (being updated)
      if (VariableFunction %in% names(private$ModuleList)){
        private$ModuleList <- private$ModuleList[private$ModuleList != VariableFunction]
        private$nodeList <- private$nodeList[private$nodeList$Calc != VariableFunction,]
      }
      aNewVariableModule <- VariableModule$new(self,VariableFunction, AggrBy = AggrBy, AggrFun = AggrFun)
      private$storeNodes(aNewVariableModule)
      private$ModuleList[[VariableFunction]] <- aNewVariableModule
      invisible(aNewVariableModule)
      
    },
    
    #' @description create and set the solver (there can only be one)
    #' @param SolverFunction The (name of) the solver defining function 
    #' @param ... passed on to init of SolverModule i.e. to CalcGraphModule$MoreParams
    #' for use in the defining function
    NewSolver = function(SolverFunction, ...){
      #existing function?
      stopifnot("function" %in% class(match.fun(SolverFunction)))
      private$solver <- SolverModule$new(self, SolverFunction, ...)
    },
    
    #' @description Run the matrix with exchange constants; 
    #' the exact calculation is described in the function that defined the solve module.
    #' This can be solving for the steady state or running the system for a period of time.
    #' Results are available 
    #' @param emissions data.frame-ish with columns "Abbr" (state defining character) and 
    #' "emis" numbers
    #' @param needdebug if T the defining function will open in debugging modus
    Solve = function(emissions, needdebug = F){
      if (is.null(private$solver)) stop ("No active solver")
      private$solver$PrepKaasM()
      if (!"data.frame" %in% class(emissions)) 
        stop("emissions should be a dataframe(like) with columns")
      if (!all(c("Abbr", "Emis") %in% names(emissions))) 
        stop("emissions should contain Abbr and Emis as columns")
      private$solver$PrepemisV(emissions)
      # the solver does the actual work
      Solution = private$solver$execute(needdebug)
      private$UpdateDL(Solution)
    },
    
    #' @description Obtain the names of the variables and tablename in which the data resides 
    #' Not needed in normal use
    #' see 
    metaData = function(){
      AllTables <- names(private$SB4Ndata)
      res <- data.frame(
        Tablenames = rep(AllTables, sapply(private$SB4Ndata, ncol)),
        AttributeNames = unname(unlist(sapply(private$SB4Ndata, names))
        ))
      res[!(res$Tablenames %in% c("Flows", "MatrixSheet")),]
    },
    
    #' @description Obtain the data for a SBvariable or a flow
    #' @param varname the name of the variable. Returns a list of variables if varname == "all"
    fetchData = function(varname="all"){
      private$FetchData(varname)
    },
    
    #' @description Pseudo-inherit data; for instance from compartment to subcompartments.
    #' Not needed in normal use.
    #' @param fromData the name of the variable before copying
    #' @param toData the name of the variable after copying
    doInherit = function(fromData, toData) {
      private$DoInherit(fromData, toData)
    },
    
    
    #' @description returns the state for a state-name. Injected from States property
    #' @param abbr the name of the state
    findState = function (abbr){
      private$States$findState(abbr)
    },
    
    #' @description returns the table, determined by the dimensions et al.
    #' @param KeyNames the dim (of The3D etc) determine which table
    whichDataTable = function(KeyNames){
      private$WhichDataTable (KeyNames)
    },
    
    #' @description runs (or tries to) the calculation graph, either the whole graph or 
    #' limited to a single process
    #' @param aProcessModule the name of the state
    #' @param mergeExisting normally leaves all other processes unchanged, but if 
    #' not mergeExisting, all process results (kaas) are cleared first
    UpdateKaas = function(aProcessModule = NULL, #(Default) NULL means calculate CalcTreeBack
                          mergeExisting = T){ 

      if (is.null(aProcessModule)) {
        NewKaas <- private$CalcTreeBack(aProcessModule = NULL)
      } else {
        if (! "ProcessModule" %in% class(aProcessModule)){ # assuming it's a string with the name of
          aProcessName <- aProcessModule
          aProcessModule <- private$ModuleList[[aProcessName]]
        }
        if (! "ProcessModule" %in% class(aProcessModule)){ # now it's an error
          stop(paste("unknown process", aProcessName))
        }
        
        NewKaas <- private$CalcTreeBack(aProcessModule)
      }

      if (is.null(private$SBkaas) | !mergeExisting){
          private$SBkaas <- NewKaas[,c("i","j","k","process")]
      } else { #merge with existing kaas; update or append if new
          Processes2Update <- unique(NewKaas$process)
          private$SBkaas <- private$SBkaas[!private$SBkaas$process %in% Processes2Update,]
          private$SBkaas <- merge(NewKaas[,c("i","j","k","process")], private$SBkaas, all = T)
      }
      #return only for purpose of transparent update; side effect is done
      invisible(NewKaas)
    },
    
    #' @description runs (or tries to) the calculation for aVariable and stores the results.
    #' After this, the SBvariable can be viewed with fetchData and this data will be used 
    #' when calculating processes, or other variables.
    #' @param aVariable the name of the Variable
    CalcVar = function(aVariable){
      TestTree <- private$nodeList[private$nodeList$ModuleType %in% c("Variable", "Flow"),]
      if (!aVariable %in% TestTree$Calc) stop(paste (aVariable, "not found"))
      private$UpdateDL(VarFunName = aVariable)
    },
    
    #' @description set a constant in the internal data, to enable use by SB variable etc.
    #' @param ... named value
    SetConst = function(...){
      #test for length == 1, name...
      private$UpdateDL(...)
    },
    
    #' @description runs (or tries to) the calculation for Variables,
    #' and continues from there to update all processes and variables 
    #' that have a depency of any of the Variables, recursively.
    #' @param Variables the name(s) of the Variable(s) that would be "dirty"
    UpdateDirty = function(Variables){#recalc DAG, starting from vector of Variables
      NotUsed <- Variables[!Variables %in% private$nodeList$Params]
      if (length(NotUsed)> 0) warning(do.call(paste, as.list(
            c("Not all Variables are used:", NotUsed))))
      private$CalcTreeForward(Variables[Variables %in% private$nodeList$Params])
    },
    
    #' @description Verifies the presence of needed variables for the calculation of 
    #' all processes
    whichUnresolved = function(){
      private$CheckTree()
    },

    #' @description removes variables appearently not needed in the calculation DAG
    #' starting from Varname upwards
    #' @param VarName the name of the Variable which inputs are cleaned, recursively
    CleanupCalcGraphAbove = function(VarName){
      private$cleanupCGAbove (VarName) 
    },
    
    #' @description Replaces a complete table in the internal data system. 
    #' Use with care.
    #' @param UpdateDF the new table, replcing the previous
    #' @param keys dims that determines the table to be updated, see whichDataTable
    #' @param TableName alternative to keys: provide the tablename. Use with care!
    UpdateData = function(UpdateDF, keys, TableName = NULL){
      if (is.null(TableName)) {
        TableName <- private$WhichDataTable(names(keys))
        stopifnot(!is.null(TableName))
      } 
      nowTable <- private$SB4Ndata[[TableName]]
      #remove current rows for keys; exception for Globals: keys = T
      if (length(keys)==1 && keys){
        private$SB4Ndata[[TableName]] <- UpdateDF
      } else {#keys of the dimensions
        DelRows <- sapply(names(keys), function(filter){
          nowTable[,filter] == keys[[filter]]
        })
        if (ncol(DelRows) > 1) {
          DelRows <- apply(DelRows, 1, all)
        }
        private$SB4Ndata[[TableName]] <- rbind(
          nowTable[!DelRows,],
          UpdateDF[,names(nowTable)])
      }
    },
    
    #' @description derive transfer processes from datalayer
    #' @param processName the one you are looking for, or (default) for "all" processes 
    FromDataAndTo = function(processName = "all"){
      #Get all form-to state combination, given
      #[3D]Processes sheets and
      #Processes columns in [3D]Sheet
      #And restrict to existing states
      PrepMatrix <- private$FetchData("Compartment")
      #do not join on SubCompart, no expansion needed there
      names(PrepMatrix) <- c("ExpandSubCompart", "compartment")
      The3Process <- lapply(The3D, function (oneOf3) {
        DTable <- private$SB4Ndata[[(paste0(oneOf3,"Processes"))]]
        if (processName != "all") {
          DTable <- DTable[DTable$process == processName,]
        }
        if (nrow(DTable) > 0) {
          #expand matrix to subcompartments?
          if (all(DTable$from %in% PrepMatrix$compartment)){
            possExpand <- lapply(1:nrow(DTable), function (xrow){
              merge(DTable[xrow,], PrepMatrix, by.x = "from", by.y = "compartment")
            }) 
            #rbind and clean columnames
            DTable <- do.call(rbind, possExpand)
            DTable$from <- DTable$ExpandSubCompart
            DTable$ExpandSubCompart <- NULL
          }
          #dito expansion for to; make a function??
          if (all(DTable$to %in% PrepMatrix$compartment)){
            possExpand <- lapply(1:nrow(DTable), function (xrow){
              merge(DTable[xrow,], PrepMatrix, by.x = "to", by.y = "compartment")
            }) 
            #rbind and clean columnames
            DTable <- do.call(rbind, possExpand)
            DTable$to <- DTable$ExpandSubCompart
            DTable$ExpandSubCompart <- NULL
          }
          
        }
        
        colnames(DTable)[colnames(DTable)=="from"] <- paste0("from",oneOf3)
        colnames(DTable)[colnames(DTable)=="to"] <- paste0("to",oneOf3)
        DTable
      })
      
      names(The3Process) <- The3D
      The3Sheet <- lapply(The3D, function (oneOf3) {
        private$SB4Ndata[[(paste0(oneOf3,"Sheet"))]]
      })
      names(The3Sheet) <- The3D
      processes3D <-lapply(The3D, function(D){
        
        #The other two dimensions; readability
        Others <- which(names(The3Sheet)!=D,arr.ind = T)
        Others1 <- Others[1]
        Others2 <- Others[2]
        AllNames1 <- lapply(The3Sheet[Others1],FUN = names)
        AllNames2 <- lapply(The3Sheet[Others2],FUN = names)
        
        #loop over processes within a dimension
        #1) list of processes
        Processes <- unique(The3Process[[D]][,"process"])
        #2 restrictions on expand?
        #search process attribute, if present indicating 'T|F' (only F matters, actually)
        # else include all elements of the dimension
        lapply(Processes, function(p){
          
          Df0 <- The3Process[[D]][The3Process[[D]]$process == p,]
          
          Oth1 <- The3Sheet[[Others1]][[The3D[Others1]]]
          if(p %in% unlist(AllNames1)){
            Oth1T <- which(The3Sheet[[Others1]][[p]] != "F")
            Oth1 <- Oth1[Oth1T]
          }
          #append the from-to
          Df1 <- data.frame(
            from = Oth1,
            to = Oth1, stringsAsFactors = F
          )
          names(Df1) <- sapply(names(Df1),paste0,The3D[Others1])
          Oth2 <- The3Sheet[[Others2]][[The3D[Others2]]]
          if(p %in% unlist(AllNames2)){
            Oth2T <- which(The3Sheet[[Others2]][[p]] != "F")
            Oth2 <- Oth2[Oth2T]
          }
          #append the from-to
          Df2 <- data.frame(
            from = Oth2,
            to = Oth2, stringsAsFactors = F
          )
          names(Df2) <- sapply(names(Df2),paste0,The3D[Others2])
          Step1Expand <- expand.grid.df(Df0, Df1, Df2)
          #remove processes with F in ProcessName - column of 2D-combination-tables
          D2combn <- combn(The3D,2)
          dummy <- apply(D2combn, 2, function(D2v){
            TableName <- paste0(do.call(paste0,as.list(D2v)),"Data")
            D2others <- private$SB4Ndata[[TableName]]
            if (p %in% colnames(D2others)) {
              ToDel <- D2others[D2others[,p]=="F" & !is.na(D2others[,p]), c(D2v[1],D2v[2])]
              for (k in 1:nrow(ToDel)){
                which1 <- which(D2others[,Others1] == ToDel[k,D2v[1]])
                which2 <- which(D2others[,Others2] == ToDel[k,D2v[2]])
                bothwhiches <- which1 & which2
                if (length(bothwhiches) < 0)
                  Step1Expand <- Step1Expand[-bothwhiches,]
              }
            }
          })
          Step1Expand
        })
      })#next D
      #flatten processes3D; 1) rbind nested level 2 2) rbind the list of dataframes
      processes3DUnnest <- lapply(processes3D, function (InList) {
        if (length(InList) == 0) return(NULL)
        if (length(InList) == 1) return(InList[[1]])
        #else rbind
        do.call(rbind,InList)
      })
      AllKaasCalc <- do.call(rbind,processes3DUnnest)
      
      #which states exist?
      exist.from <- apply(AllKaasCalc, 1, function(othrow) {
        any(private$States$asDataFrame$Scale == othrow["fromScale"] &
              private$States$asDataFrame$SubCompart == othrow["fromSubCompart"] &
              private$States$asDataFrame$Species == othrow["fromSpecies"])
      })
      exist.to <- apply(AllKaasCalc, 1, function(othrow) {
        any(private$States$asDataFrame$Scale == othrow["toScale"] &
              private$States$asDataFrame$SubCompart == othrow["toSubCompart"] &
              private$States$asDataFrame$Species == othrow["toSpecies"])
      })
      return(AllKaasCalc[exist.from & exist.to, ]) 
      
    }
    
  ), 
  
  active = list(
    #' @field states getter for r.o. property
    states = function(value) {
      if (missing(value)) {
        private$States
      } else {
        stop("`$states` are set by new()", call. = FALSE)
      }
    },
    #' @field kaas getter for r.o. property (all k's)
    kaas = function(value) {
      if (missing(value)) {
        #create join met states
        if (is.null(private$SBkaas)) return(NULL)
        #else
        private$ijAddState(private$SBkaas)
      } else {
        warning("`$kaas` should be set by UpdateKaas(), format is strict!", call. = FALSE)
        private$SBkaas <- value
      }
    },
    #' @field nodelist getter for r.o. property
    nodelist = function(value) {
      if (missing(value)) {
        private$nodeList
      } else {
        stop("use the $NewProcess(), $NewCalcVariable() and $NewFlow() methods to construct a node-list", call. = FALSE)
      }
    },
    #' @field moduleList getter for r.o. property
    moduleList = function(value) {
      if (missing(value)) {
        private$ModuleList
      } else {
        stop("use the $NewProcess(), $NewCalcVariable() and $NewFlow() methods to construct a node-list", call. = FALSE)
      }
    }
  ),
  
  private = list(
    SB4Ndata = NULL,
    States = NULL,
    substance = NULL,
    SBkaas = NULL,
    ModuleList =  NULL,
    nodeList = NULL,
    solver = NULL,
    l_postPoneList = NULL,

    DoInherit = function(fromDataName, toDataName){
      #Inherits from global, Matrix, compartment, or a subset of dimensions The3D of the toData
      # parameters should be a fetch-able string
      fromData <- private$FetchData(fromDataName)
      
      toData <- self$fetchData(toDataName)
      stopifnot("data.frame" %in% class(toData))
      #toData should have 1 - 3 dimension
      toDataDims <- The3D %in% names(toData)
      stopifnot(any(toDataDims), sum(toDataDims) <= 3)
      fromTableNames <- c(The3D, "Matrix", "Compartment")
      fromDataDims <- fromTableNames %in% names(fromData)
      names(fromDataDims) <- fromTableNames
      # Include the NA rows (removed by fetchData)
      DimensionNames <- The3D[The3D %in% names(toData)]
      #The name
      FullScaffold <- self$whichDataTable(KeyNames = DimensionNames)
      #The actual data.frame
      FullScaffold <- private$SB4Ndata[[FullScaffold]][,DimensionNames, F] #don't drop frame
      #extend toData to the Full
      toData <- left_join(FullScaffold, toData)

      if (!any(fromDataDims)){
        #it's a constant, we hope?
        stopifnot ((is.numeric(fromData) && length(fromData) == 1) ) 
        toData$toRename = fromData
        names(toData)[names(toData)=="toRename"] <- fromDataName
        #toData
      } else { #it should be a data.frame
        stopifnot("data.frame" %in% class(fromData))
        if (!fromDataDims["Matrix"]) {
          toData <- left_join(toData, fromData)
        } # else #The from is already inherited to Subcompartment at initialisation
      }
      updateRows <- is.na(toData[,toDataName])
      toData[updateRows, toDataName] <- toData[updateRows, fromDataName]
      NeededNames <- names(toData)[names(toData) %in% c(The3D, toDataName)] 
      toData <- toData[!is.na(toData[,toDataName]), NeededNames]
      private$UpdateDL(VarFunName = toData)
      
    },
    
    storeNodes = function(aNewModule){
      moduleType <- switch (class(aNewModule)[1],
        "VariableModule" = "Variable",
        "FlowModule" = "Flow",
        "ProcessModule" = "Process"
      )
      aNeedVars <- unique(aNewModule$needVars) #there might be a from. AND a to. parameter!
      if (length(aNeedVars) == 0){
        NewNodes <- data.frame(Calc = aNewModule$myName, Params = "", ModuleType = moduleType)
      } else NewNodes <- data.frame(Calc = aNewModule$myName, Params = aNeedVars, ModuleType = moduleType)
      if (is.null(private$nodeList)) {
        private$nodeList <- NewNodes
      } else private$nodeList <- rbind(private$nodeList, NewNodes)

    },
    
    cleanupCGAbove = function(VarName) {
      browser()
      VarNameVector <- VarName
      MoreVarNames <- VarName
      repeat {
        candidotes <- private$nodeList[private$nodeList$Calc %in% MoreVarNames,]
        NoCandid <- private$nodeList[private$nodeList$Params %in% candidotes$Params &
                                       !private$nodeList$Calc %in% VarNameVector,]
        MoreVarNames <- candidotes$Params[!candidotes$Params %in% NoCandid$Params]
        if (length(MoreVarNames) == 0) break
        VarNameVector <- unique(c(VarNameVector, MoreVarNames))
      }
      #TODO remove from data?
      cat (do.call(paste, c("removed nodes for", VarNameVector)))
      private$nodeList <- private$nodeList[private$nodeList$Calc %in% VarNameVector,]
    },
    
    ijAddState = function(ijTable) {
      #local and specific function to rename WasNames, adding PrePart if appropp 
      PreNames <- function(WasNames, PrePart) {
        for(WNi in 1:length(WasNames)){ #this doesn't take time
          if (WasNames[WNi] %in% c("Scale","SubCompart","Species","Abbr" ))
            WasNames[WNi] <- paste(PrePart, WasNames[WNi], sep = "")
        }
        return(WasNames)
      }
      fromS <- private$States$asDataFrame[ijTable$i,]
      names(fromS) <- PreNames(names(fromS), "from")
      toS <- private$States$asDataFrame[ijTable$j,]
      names(toS) <-  PreNames(names(toS), "to")
      #only 1 of 3D is different, leave it out for clarity?
      toS$toSpecies [toS$toSpecies == fromS$fromSpecies] <- ""
      toS$toSubCompart [toS$toSubCompart == fromS$fromSubCompart] <- ""
      toS$toScale [toS$toScale == fromS$fromScale] <- ""
      return(cbind(ijTable, fromS[,-match("i", names(fromS))], toS[,-match("i", names(toS))]))
    },
    
    CheckTree = function(){ #recursive checking the calculation of the DAG of processes/variables
      TestTree <- private$nodeList
      #remove the ones present in the data; these you have
      MetaData <- self$metaData()
      TestTree$Params[TestTree$Params %in% MetaData$AttributeNames] <- ""
      AllWant <- unique(TestTree$Calc)
      TestTree$Check <- TestTree$Params == "" 
      repeat{
        CountCantdo <- length(which(!TestTree$Check))
        CantDo <- unique(TestTree$Calc[!TestTree$Check])
        CanDo <- AllWant[!AllWant %in% CantDo]
        TestTree$Check <- TestTree$Check | TestTree$Params %in% CanDo
        if (CountCantdo == length(which(!TestTree$Check))) break
      }
      unique(TestTree$Params[!TestTree$Check])
    },
    
    CalcTreeForward = function(DirtyVariables){ #calculation of variables and kaas
      if (is.null(DirtyVariables)) stop("Cannot CalcTreeForward without a(starting/dirty)Variable")
      TestTree <- private$nodeList
      browser()
      AllCan <- TestTree$x
    },
    
    CalcTreeBack = function(aProcessModule){ #calculation of variables and kaas
     
      #treat the objects that do not use private$ModuleList separate
      if ("ClassicNanoProcess" %in% class(aProcessModule))   {
        return(aProcessModule$execute())
        
      } else {
        kaaslist <- list() # all the dataframes of (future) new kaas
        if (is.null(aProcessModule)) {
          #exception for the CalcGraph* that are calculated after the tree i.e. molecular deposition special
          TestTree <- private$nodeList[!private$nodeList %in% unlist(private$l_postPoneList)]
        } else {
          if (! "ProcessModule" %in% class(aProcessModule)){ # assuming it's a string with the name of
            aProcessModule <- private$nodeList[[aProcessModule]]
          }
          TestTree <- private$nodeList[private$nodeList$ModuleType %in%  c("Variable", "Flow") | 
                                         private$nodeList$Calc==aProcessModule$myName,]
        }
        AllWant <- unique(TestTree$Calc)
        MetaData <- self$metaData()
        TestTree$Params[TestTree$Params %in% MetaData$AttributeNames] <- ""
        
        repeat{
          TestTree$Check <- TestTree$Params == ""
          CountCantdo <- length(which(!TestTree$Check))
          CantDo <- unique(TestTree$Calc[!TestTree$Check])
          CanDo <- AllWant[!AllWant %in% CantDo]
          if (length(CanDo > 0)){
            for (i in 1:length(CanDo)){
              CalcMod <- private$ModuleList[[CanDo[i]]]
              if ("VariableModule" %in% class(CalcMod) | "FlowModule" %in% class(CalcMod)) { #update DL
                succes <- private$UpdateDL(CanDo[i])
                if (nrow(succes) < 1) warning(paste(CanDo[i],"; no rows calculated"))
              } else { # a process; add kaas to the list
                
                kaaslist[[CalcMod$myName]] <- CalcMod$execute()
              }
              TestTree$Params[TestTree$Params == CanDo[i]] <- ""
              TestTree <- TestTree[-which(TestTree$Calc == CanDo[i]),] 
              AllWant <- AllWant[-which(AllWant == CanDo[i])]
            }
          } else {
            #anything left; enough?
            if (length(CantDo)>0){
              if(is.null(aProcessModule)){
                CantDoProcesses <- CantDo[private$ModuleList[CantDo]$ModuleType == "Process"]
                stop(do.call(paste(c(list("Can't calculate"), as.list(CantDoProcesses)))))
              } else {
                #MissingVar <- TODO
                stop(paste("Can't calculate", aProcessModule$myName))
              }
            }
            break #repeat
          }
        }
        # do postponed
        for (postNames in c(private$l_postPoneList["VarFunctions", "FlowFunctions", "ProcesFunctions"])){ #force that order
          CalcMod <- private$ModuleList[[CanDo[i]]]
          if ("VariableModule" %in% class(CalcMod) | "FlowModule" %in% class(CalcMod)) { #update DL
            succes <- private$UpdateDL(CanDo[i])
            if (nrow(succes) < 1) warning(paste(CanDo[i],"; no rows calculated"))
          } else { # a process; add kaas to the list
            kaaslist[[CalcMod$myName]] <- CalcMod$execute()
          }
        }
      }
      #select kaas names only for all kaaslist - data.frames
      kaaslist <- lapply(kaaslist, 
                         dplyr::select, c("process", "fromScale","fromSubCompart","fromSpecies",
                                      "toScale","toSubCompart","toSpecies","k"))
      VolleKaas <- do.call(rbind, kaaslist)
      
      if(length(VolleKaas) > 0) { #length is 0 if no processes 
        #format to old-school k names; merge with existing kaas and return a diff
        # 1 match from and to states: i and j 
        return( data.frame(
          i = private$FindStatefrom3D(data.frame(
            Scale = VolleKaas$fromScale,
            SubCompart = VolleKaas$fromSubCompart,
            Species = VolleKaas$fromSpecies
          )),
          j = private$FindStatefrom3D(data.frame(
            Scale = VolleKaas$toScale,
            SubCompart = VolleKaas$toSubCompart,
            Species = VolleKaas$toSpecies
          )),
          k = VolleKaas$k,
          process = VolleKaas$process
        ) )
      } else {return(NULL)}
      
    },

    #find the indices + states from a dataFrame with columns "Scale" "SubCompart" "Species"   
    FindStatefrom3D = function(df3Ds){
      stopifnot(all(The3D %in% names(df3Ds)))
      ret <- sapply(1:nrow(df3Ds), function(i){
        which(private$States$asDataFrame$Scale == df3Ds$Scale[i] &
                private$States$asDataFrame$SubCompart == df3Ds$SubCompart[i] &
                private$States$asDataFrame$Species == df3Ds$Species[i])
      })
      if (anyNA(ret)) {
        cat(df3Ds[is.na(ret), c("Scale", "SubCompart", "Species")])
        stop("State(s) not found")
      } else {
        return(ret)
      }
    },
    
    #facilitate access to self$SB4Ndata
    
    #' description fetch the values for a dataframe with parameters/dimensions (any purpose),
    #' param varname name of variable to find
    ## '  @param ParamName.Dimensions data.frame with columns ParamName and the 3Dimensions() use merge x.all = T
    #' return vector values (if found; else NA)
    FetchData = function(varname) {
      #hack; because of ugly deposition dependency
      if (varname == "kaas") {
        return(self$kaas)
      }
      #Now for the actual fetching of variables; either calculated or data
      if (varname %in% names(private$SB4Ndata)) {
        #give the full table
        return(private$SB4Ndata[[varname]])
      } else {
        #is it a flow?
        fromFlows <- private$SB4Ndata[["Flows"]]
        if (any(fromFlows$FlowName == varname)) {
          fromFlows <- fromFlows[fromFlows$FlowName == varname,]
          attr(fromFlows, which = "isflow") <- T
          return(fromFlows)
        } else {
          QSARnames <- names(private$SB4Ndata[["QSARtable"]])
          if (varname %in% QSARnames) {
            ChemClassNow <- private$SB4Ndata[["Globals"]]$ChemClass
            res = private$SB4Ndata[["QSARtable"]][private$SB4Ndata[["QSARtable"]]$QSAR.ChemClass == ChemClassNow,varname]
            return(varname = res)
          } else {
            #normal behaviour: find the attribute and return its "view"
            MetaData <- self$metaData()
            if (varname=="all"){ #return an overview (listing of variables)
              return( #all flows and all attribute names excluding key fields
                sort(c(
                  unique(fromFlows$FlowName),
                  QSARnames,
                  MetaData$AttributeNames[!MetaData$AttributeNames %in% c(
                    The3D, paste("to.", The3D, sep = ""),
                    "Default", "from", "to", "process", "AbbrS"
                  ) & ! MetaData$Tablenames %in% c("SubstanceCompartments")]
                ))
              )
            }
          }
          
          Attrn <- MetaData$Tablenames[which(varname == MetaData$AttributeNames)]
          
          if (length(Attrn)==0) {
            allVars <- MetaData[!MetaData$AttributeNames %in% c(
              The3D, paste("to.", The3D, sep = ""),
              "Default", "from", "to", "process", "AbbrS"
            ) & ! MetaData$Tablenames %in% c("SubstanceCompartments"),]
            grepVars <- do.call(paste,
                                as.list(allVars$AttributeNames[grep(varname, allVars$AttributeNames, ignore.case=TRUE)]))
            stop (paste("Cannot find property", varname, "; but found", grepVars))
          } 
          
          #Attrn <- levels(Attrn)[Attrn] R Version 4 now 
          if(length(Attrn) > 1) stop (paste0(varname, "in multiple tables:", Attrn))
          
          subvec <- MetaData[MetaData$Tablenames == Attrn, "AttributeNames"]
          Dims <- c(The3D, paste("to.", The3D, sep = "")) %in% subvec
          names(Dims)[Dims] <- c(The3D, paste("to.", The3D, sep = ""))[Dims]
          DimsVarCols <- c(names(Dims)[Dims], varname)
          
          #Check if unit should be converted to SI
          unitTable <- private$SB4Ndata[["Units"]]
          SIexpression <- unitTable$ToSI[unitTable$VarName == varname]
          
          isExpression <-  length(SIexpression) == 1 &&
            !is.na(SIexpression) && 
            SIexpression != "" && 
            SIexpression != "NA"
          
          Doexpression <- function (varname, x, SIexpression){ #execute expression to convert to SI
            #NB varname is local here
            assign(varname, x)
            eval(parse(text = SIexpression))
          }
          if (any(Dims)) { #it's a data.frame with varname as column
            res <- private$SB4Ndata[[Attrn]][,DimsVarCols]
            #remove the NA's from the pivot_wider
            ToDel <- is.na(res[,varname])
            
            if (isExpression){
              # ADJUST TO SI UNITS!!!
              
              res[,varname] <- Doexpression(varname, res[,varname], SIexpression)
            }
            return(res[!ToDel, DimsVarCols])
          } else { #it's atomic 
            res <- private$SB4Ndata[[Attrn]][,DimsVarCols]
            if (isExpression){
              #  ADJUST TO SI UNITS!!!
              res <- Doexpression(varname, res, SIexpression)
            }
            return(varname = res)
          }
        }
      }
    },
    
    #' #description calculate and update a variable in the "datalayer" (SB4Ndata-dataframes)
    #' #param TableNm the name of the data.frame in SB4Ndata
    #' #param VarFunName the name of the variable to be updated AND the name of the function
    #' #param DIMRestrict optional list of restrictions; each element is c(DimName, Comparator, Value)
    #' #return side-effect; vector?
    UpdateDL = function(VarFunName = NULL, DIMRestrict = NULL, ...) {
     
      MetaData <- self$metaData()
      if (is.null(VarFunName)) {
        inp <- list(...)
        VarFunName <- names(inp)
        NewData <- inp
        Attrn <- MetaData$Tablenames[match(VarFunName, MetaData$AttributeNames)]
      } else {
        if ("data.frame" %in% class(VarFunName)) { 
          NewData <- VarFunName
          VarFunName <- names(NewData)[length(NewData)]
          Attrn <- MetaData$Tablenames[match(VarFunName, MetaData$AttributeNames)]
        } else { #
          Attrn <- MetaData$Tablenames[match(VarFunName, MetaData$AttributeNames)]
          #exist VarFunName as VariableModule?
          if (! VarFunName %in% names(private$ModuleList)) stop(paste("Can't find", VarFunName, "as VariableModule"),call. = F)
          VarFun <- private$ModuleList[[VarFunName]]
          NewData <- VarFun$execute()
          #adjust column names 
          if ("FlowModule" %in% class(VarFun)  ) {
            #to conform the full Flows table)
            NewData$FlowName <- VarFunName
            if (!"toScale" %in% names(NewData)){
              NewData$toScale <- NewData$fromScale
            }
            if (!"toSubCompart" %in% names(NewData)){
              NewData$toSubCompart <- NewData$fromSubCompart
            }
          } 
          if ("VariableModule" %in% class(VarFun)){
            if(is.null(names(NewData))){
              names(NewData) <- VarFunName
            }
          }
        }
      }
      
      # where should it go? flows or Where the dimensions are
      # and prep diff 
      if ("flow" %in% names(NewData)){
        Target.Table <- "Flows"
        if (VarFunName %in% private$SB4Ndata[["Flows"]]$FlowName){
          diffTable <- private$FetchData(VarFunName)
          names(diffTable)[names(diffTable)==VarFunName] <- paste("old",VarFunName,sep = "_")
        }
      } else {
        Target.Table <- private$WhichDataTable(names(NewData))
        if (!is.na(Attrn)){ #is already in the DL
          diffTable <- private$FetchData(VarFunName)
          names(diffTable)[names(diffTable)==VarFunName] <- paste("old",VarFunName,sep = "_")
        }
      }
      #Special Case a "constant" or substance property
      if (Target.Table == "Globals") {
        private$SB4Ndata[[Target.Table]][names(NewData)] <- NewData
      } else {
        #delete and merge the DL-table
        if (VarFunName %in% names(private$SB4Ndata[[Target.Table]])){
          numCol <- match(VarFunName, names(private$SB4Ndata[[Target.Table]]))
          private$SB4Ndata[[Target.Table]] <- private$SB4Ndata[[Target.Table]][,-numCol]
        } 
        private$SB4Ndata[[Target.Table]] <- merge(private$SB4Ndata[[Target.Table]], NewData, all = T)
      }
      
      #just to show, side effect is in DL
      if (!exists("diffTable")){ #new in the DL
        return(NewData)
      } else {
        return(merge(diffTable, NewData, all = T)) 
      }
    },
    
    WhichDataTable = function(KeyNames){# Which table ("1D"sheet or "2-3D"Data) in SB4N.data has identical dimensions?
      # Special case, a Flows has special keynames for easy merging when calculation processes; 
      # none of the actual keynames are in The3D
      if (length(KeyNames) == 0){ #either Globals or a flow
        if ("flow" %in% KeyNames){
          return("Flows")
        }
      }
        
      MetaData <- self$metaData()
      MetaDims <- data.frame(
        SB4N.Table = unique(MetaData$tablenm)
      )
      WhichD <- The3D[The3D %in% KeyNames]
      AppendString <- switch (length(WhichD)+1,
        "Globals",
        "Sheet",
        "Data",
        "Data")
      do.call(paste,as.list(c(WhichD, AppendString, sep = "")))
      
    }

  )
)

#' #S3 function TODO print method? see docs R6
#' #' @export
#' Summary.SBcore <- function(MySBcore) {
#'   stopifnot(all(sapply(c("kaas", "states"), `%in%`, names(SBWorld))))
#'   Scal <- table(SBWorld$states$Scale)
#'   Scales  <- Scal[Scal > 0]
#'   SubC <- table(SBWorld$states$SubCompart)
#'   SubCompartments <- SubC[SubC > 0]
#'   Specy <- table(SBWorld$states$Species)
#'   Species <- Specy[Specy > 0]
#'   BoxDimensions <- list(Scales = Scales, 
#'                         SubCompartments = SubCompartments,
#'                         Species = Species)
#'   processM <- GetProcesses(SBWorld$states)
#'   processkaas <- aggregate(fromSubCompart~process, 
#'                            data = processM, FUN = length)
#'   names(processkaas) <- c("process", "No.Kaas" )
#'   summlog10kaas <- table(round(log10(SBWorld$kaas$k)))
#'   list(BoxDimensions = BoxDimensions,
#'        processKaas = processkaas, log10kaasTafel = summlog10kaas)
#'   
#' }