#' @title SBcore
#' @description the centre object for running simple box with an R script.
#' @import R6
#' @export
SBcore <- R6::R6Class("SBcore",
  public = list(
    #' @description init
    #' @param SBooDataLocation The path to data for SB
    #' @param SuBmode Either string ("Molecular"|"Particle") or the filepath to 
    #' (substance)data for SB (overruling any default data).
    #'  SBmode will be set to "Molecular"|"Particle" accordingly
    initialize = function(SBooDataLocation, SuBmode = "Molecular", debugR = FALSE){
      
      private$myDataLocation <- SBooDataLocation
      private$debugR <- debugR
      self$SuBmode <- SuBmode # also sets private$SBmode

      self$RawTables <- ReadRawData$new(private$myDataLocation)
      
      Rules <- read.csv(file.path(SBooDataLocation, "StateRules.csv"))
      ModeRules <- read.csv(file.path(SBooDataLocation, "StateModeRules.csv")) 
      
      private$States <- private$getStates(
        Scales = self$RawTables$DimsTables[["Scale"]],
        Species = self$RawTables$DimsTables[["Species"]],
        SubComparts = self$RawTables$DimsTables[["SubCompart"]],
        Rules = Rules,
        ModeRules = ModeRules)
      
      private$States$myCore <- self

      # init reactives
      private$myReactiveDAG <- ReactiveDAG$new(debugR = private$debugR, parent = self)

      # and fill it with the data from Defs
      private$data2DAG(self$RawTables$CONST)
      for (DefTable in self$RawTables$VarsTables){
        private$data2DAG(DefTable, Table_type = "Vars")
      }
      for (DefTable in self$RawTables$ToWidenTables){
        private$data2DAG(DefTable, Table_type = "ToWiden")
      }
      
      if (is.list(self$SuBmode)) {
        # set substance properties from list/json, all but SBmode
        SubstanceList <- self$SuBmode$Substance[!(names(self$SuBmode$Substance) == "SBmode")]
        private$Block2DAG(SubstanceList)
        NonSubstanceList <- self$SuBmode[!(names(self$SuBmode) == "Substance")]
        private$Block2DAG(SubstanceList)
      }

      # used as "setting" in ReactiveDAG
      speciesInUse <- private$States$asDataFrame |>
        dplyr::pull(Species) |> unique()
      private$myReactiveDAG$add_source("species", speciesInUse)
      
      # flows 4 advection
      flows <- self$RawTables$FlowIO
        
      # (other) processes
      process_functions = read.csv(file.path(SBooDataLocation, "Processes4SpeciesTp.csv"))
      processes4Mode <- process_functions$Process[process_functions[, private$SBmode] == "X"]

      SBvars <- self$build_DAG(unique(flows$FlowName), processes4Mode)
      
      message("processes (besides advection) ", do.call(paste, as.list(processes4Mode)))
      
      # assuming unknown variable should be NA (function should raise the error)
      SBvars$missing_vars
      misvarstr <- do.call(paste, SBvars$missing_vars)
      lapply(SBvars$missing_vars, function(x) private$myReactiveDAG$add_source(x, NA))
      message(paste("unknown parameters (set to NA): ", misvarstr))
      
      # store (names of) flows and processes for UpdateKaas()
      private$flows <- unique(flows$FlowName)
      # warning for updated Advection?
      if ("Advection" %in% processes4Mode) {
        warning("no need to define Advection")
        processes4Mode <- processes4Mode[processes4Mode != "Advection"]
      } 
      private$processes <- paste0("k_", processes4Mode)
      
    },

    #' @description create and set the solver (there can only be one)
    #' @param SolverFunction The (name of) the solver defining function 
    #' @param ... passed on to init of SolverModule i.e. to CalcGraphModule$MoreParams
    #' for use in the defining function
    NewSolver = function(SolverFunction, ...){
      #existing function?
      stopifnot("function" %in% class(match.fun(SolverFunction)))
      private$solvername <- SolverFunction
      private$solver <- SolverModule$new(self, SolverFunction, ...)
    },
    
    set_substance = function(pathOrList) {
      
      if (!isalist(pathOrList)){
        # find file as full path, in myDataLocation or myDataLocation/substances
        if (!file.exists(pathOrList)) {
          path <- file.path(private$myDataLocation, pathOrList)
        } else {
          if (!file.exists(path)){
            path <- file.path(private$myDataLocation, "substances", pathOrList)
          } else {
            path <- pathOrList
          }
        }
        if (!file.exists(path)) {
          stop(paste("substance json file does not exist:", path))
        }
        json <- read_json(path, simplifyVector = FALSE)
        substances <- json$substances
        
        out <- list()
        
        for (sub in substances) {
          sub_id   <- sub$id   %||% NA_character_
          sub_name <- sub$name %||% NA_character_
          props    <- sub$properties
          if (is.null(props)) next
          
          for (prop_name in names(props)) {
            prop <- props[[prop_name]]
            
            # "table" if it has a 'values' array of records
            if (!is.null(prop$values) && is.list(prop$values)) {
              unit <- prop$unit %||% NA_character_
              
              df <- map_dfr(prop$values, function(row) {
                row <- row %||% list()
                tibble(!!!row)
              })
              
              df <- df |>
                mutate(
                  substance_id   = sub_id,
                  substance_name = sub_name,
                  property       = prop_name,
                  unit           = unit,
                  .before = 1
                )
              
              key <- paste(sub_id %||% sub_name %||% "unknown", prop_name, sep = "_")
              out[[key]] <- df
            }
          }
        }
        #create / update reactives
        
      }
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
    #' @description Run the matrix with exchange constants; 
    #' the exact calculation is described in the function that defined the solve module.
    #' This can be solving for the steady state or running the system for a period of time.
    #' Results are available 
    #' @param emissions data.frame-ish with columns "Abbr" (state defining character) and 
    #' "emis" numbers
    #' @param needdebug if T the defining function will open in debugging modus
    Solve = function(emissions, needdebug = F, ...){
      
      if (is.null(private$solver)) {
        warning("No active solver")
        return(NULL)
      }
      #browser()
      # prepare the kaas
      private$solver$PrepKaasM()
      
      private$solver$execute(needdebug = needdebug, emissions, ...)
      
    },
    
    #' @description Creates inverse functions for each row of the given dataframe
    #' @param paramdf A dataframe with one row for eacht distribution. Must at 
    #' least contain the columns "Distribution" (which can contain "triangular", "uniform" or "normal"), 
    #' "a", "b", and "c".
    #' 
    #' "a", "b" and "c" contain the following:
    #' triangular: a = Minimum, b = Maximum, c = Peak value
    #' normal: a = Mean, b = Standard deviation 
    #' uniform: a = Minimum , b = Maximum
    makeInvFuns = function(paramdf){
      
      # Define functions for each row based on the distribution type
      varFuns <- apply(paramdf, 1, function(aRow) {
        dist_type <- aRow["Distribution"]
        
        if (dist_type == "triangular" || dist_type == "Triangular") {
          prepArgs <- as.list(as.numeric(aRow[c("a", "b", "c")]))
          names(prepArgs) <- c("a", "b", "c")
        } else if (dist_type == "normal" || dist_type == "Normal") {
          prepArgs <- as.list(as.numeric(aRow[c("a", "b")]))
          names(prepArgs) <- c("a", "b")
        } else if (dist_type == "uniform" || dist_type == "Uniform") {
          prepArgs <- as.list(as.numeric(aRow[c("a", "b")]))
          names(prepArgs) <- c("a", "b")
        } else if (dist_type == "log uniform" || dist_type == "Log uniform") {
          prepArgs <- as.list(as.numeric(aRow[c("a", "b")]))
          names(prepArgs) <- c("a", "b")
        } else if (dist_type == "TRWP_size") {
          prepArgs <- as.list(as.character(aRow[c("d")]))
          names(prepArgs) <- c("d")
        } else if (dist_type == "weibull" || dist_type == "Weibull") {
          prepArgs <- as.list(as.numeric(aRow[c("a", "b", "c")]))
          names(prepArgs) <- c("a", "b", "c")
        } else if (dist_type == "log normal" || dist_type == "Log normal") {
          prepArgs <- as.list(as.numeric(aRow[c("a", "b", "c")]))
          names(prepArgs) <- c("a", "b", "c")
        } else if (dist_type == "power law" || dist_type == "Power law") {
          prepArgs <- as.list(as.numeric(aRow[c("a", "b", "c")]))
          names(prepArgs) <- c("a", "b", "c")
        } else if (dist_type == "trapezoidal" || dist_type == "Trapezoidal") {
          prepArgs <- as.list(as.numeric(aRow[c("a", "b", "c", "d")]))
          names(prepArgs) <- c("a", "b", "c", "d")
        } else {
          stop("Unsupported distribution type")
        }
        
        # Create the inverse CDF function using the prepared arguments
        Make_inv_unif01(fun_type = dist_type, pars = prepArgs)
      })
      
      return(varFuns)
    },
    
    #' @description Export the matrix of speed constants, aka Engine, to an excel file
    exportEngine = function(excelFile) {
      if (is.null(private$solver)) {
        warning("No active solver")
        return(NULL)
      }
      ToExport <- private$solver$PrepKaasM()
      dframe2excel(as.data.frame(ToExport), outxlsx = excelFile)
    },
    
    #'@description Export the matrix to World 
    exportEngineR = function(){
      if (is.null(private$solver)) {
        warning("No active solver")
        return(NULL)
      }
      private$solver$PrepKaasM()
    },
    
    #'@description Save the last calculated masses in the core
    K_matrix = function(){
      #browser()
      if (is.null(private$solver)) {
        stop("No active solver")
      }
      private$solver$GetK_matrix()
    },
    
    #'@description Save the last calculated masses in the core
    Masses = function(){
      #browser()
      if (is.null(private$solver)) {
        stop("No active solver")
      }
      private$solver$GetMasses()
    },
    
    VariableValues = function(){
      if (is.null(private$solver)) {
        stop("No active solver")
      }
      private$solver$GetVarValues()
    },
    
    #'@description Function to obtain steady state concentrations, using the solution saved in world.
    Concentration = function(){
      if (is.null(private$solver)) {
        stop("No active solver, (then Solve ..., then ask again)")
        return(NULL)
      }
      private$solver$GetConcentrations()
    },
    
    #'@description Save the last used emissions in the core
    Emissions = function(){
      #browser()
      if (is.null(private$solver)) {
        stop("No active solver")
      }
      private$solver$GetEmissions()
    },
    
    #' @description Return the appropriate concentration plot
    PlotConcentration = function(scale = NULL, subcompart = NULL){
      if (is.null(private$solver)) {
        stop("No active solver")
      }
      private$solver$GetConcentrationPlot(scale, subcompart)
    },
    
    #' @description Return the appropriate solution plot
    PlotMasses = function(scale = NULL, subcompart = NULL){
      if (is.null(private$solver)) {
        stop("No active solver")
      }
      private$solver$GetMassesPlot(scale, subcompart)
    },
    
    #' @description Return the appropriate mass distribution plot
    PlotMassDistribution = function(scale = NULL){
      if (is.null(private$solver)) {
        stop("No active solver")
      }
      private$solver$GetMassDist(scale)
    },
    
    #' @description Injection from SolverModule
    MassesAsRelational = function(...){
      if (is.null(private$solver)) {
        warning("No active solver")
        return(NULL)
      }
      private$solver$MassesAsRelational(...)
    },
    
    #' @description Obtain the names of the variables and tablename in which the data resides 
    #' Not needed in normal use
    #' see 
    metaData = function(){
      list(sources = names(private$myReactiveDAG$sources),
           nodes = names(private$myReactiveDAG$nodes))
    },
    
    #' @description Obtain the data for a SBvariable or a flow
    #' @param varname the name of the variable. Returns a list of variables if varname == "all"
    fetchData = function(varname="all"){
      if (varname == "all") {
        private$myReactiveDAG$knowns
      }  else {
        private$myReactiveDAG$get_value(varname)
      }
    },
    
    reactiveData = function(varname = "all"){
      if (varname == "all") {
        private$myReactiveDAG$dataNames
      }  else {
        if (varname %in% private$myReactiveDAG$dataNames){
          private$myReactiveDAG$sources[[varname]]
        } else {
          message(paste("unknown:", varname)) 
        }
      }
    },
    
    #' @description fetch specific values from core
    #' @param aVar data.frame-ish as used 
    fetch_current = function(withoutValues) {
      private$myReactiveDAG$trigger(withoutValues)
    },
    
    #' @description function to obtain the data for a variable or flow, including the units whenever present in the Units csv
    #' @param varname name of the variable
    fetchDataUnits = function(varname="all"){
      #browser()
      fd <- private$FetchData(varname)
      
      if(varname != "all"){
        unitTable <- private$SB4Ndata[["Units"]]
        if(varname %in% unitTable$VarName){
          
          # Check if 'fd' is a data.frame or an atomic vector
          if (is.data.frame(fd)) {
            unit <- unitTable |> filter(VarName == varname)
            unit <- unit$Unit
            fd <- fd |> mutate(Unit = unit)
            
          } else if (is.atomic(fd) && length(fd) == 1) {
            # If fd is an atomic vector (e.g., a named number), attach unit as attribute
            unit <- unitTable |> filter(VarName == varname)
            unit <- unit$Unit
            
            # Add the unit attribute to the vector
            attr(fd, "unit") <- unit
          }
        }
      }
      return(fd)
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
    
    #' @description limits the states according to the filter
    #' @param ... the name of one of the dimensions of the states (ScaleName, SubCompartName, SpeciesName)
    filterStatesFrame = function(inDataFrame){
      if (!"data.frame" %in% class(inDataFrame)){
        #browser()
        stop("expected a data.frame as parameter, to apply the filterstatus PROPERTY")
      }
      theFilter <- private$filterstates
      outframe <- inDataFrame
      for (argName in names(theFilter)) {
        DimIdNames <- private$FetchData(argName)
        DimName <- names(DimIdNames)[names(DimIdNames) %in% The3D]
        if (any(endsWith(names(outframe), DimName))) {
          FilterDim <- DimIdNames[DimIdNames[,argName] == theFilter[[argName]],DimName]
          likeDimNames <- names(outframe)[endsWith(names(outframe), DimName)]
          for (likeDimName in likeDimNames){
            outframe <- outframe[(!is.na(outframe[,likeDimName])) & outframe[,likeDimName] %in% FilterDim,]
          }
        }
      }
      return(outframe)
    },
    
    #' @description returns the table, determined by the dimensions et al.
    #' @param KeyNames the dim (of The3D etc) determine which table
    whichDataTable = function(KeyNames){
      private$WhichDataTable (KeyNames)
    },
    
    #' @description runs (or tries to) the calculation DAG
    #' Parameters are no longer usefull with R6ReactiveDAG
    UpdateKaas = function(){ 

      # Calculate (fetch) flows, assign to Advection process and add to kaas
      fromVolume <- self$fetchData("Volume")
      asProcess_list <- lapply(private$flows, function(fl) {
        flow_calculated <- self$fetchData(fl)
        k2_Advection(flow_calculated, fromVolume)
      })
      # Calculate (fetch) all other processes
      ProcessKaas <- lapply(private$processes, function(procs) {
        self$fetchData(procs)
      })
      
      # rename all to k / expand 2 all from/to dim
      ProcessKaas <- private$PrepareKaaslist(ProcessKaas)
      private$SBkaas <- do.call(rbind, c(ProcessKaas, asProcess_list))
      
      #return only for purpose of transparent update; side effect is done
      invisible(private$SBkaas)
    },
    
    #' @description Tries to create all non-data reactives that current processes and all flows need
    build_DAG = function(flowNames, ModeProcesses) {
      ns <- getNamespace("sboo")
      my_ls <- ls(ns)
      
      processnames <- paste("k", ModeProcesses, sep = "_")
      # exclude "k_Advection", it's now automatically added from flows 
      processnames <- processnames[processnames != "k_Advection"]
      
      stopifnot(all(processnames %in% my_ls))
      
      stopifnot(all(flowNames %in% my_ls))
      
      molecular_deposition_proces = "k_Deposition"
      
      if (molecular_deposition_proces %in% processnames){
        # obtain nasty deposition, depending on other air sinks (param "OtherkAir")
        fromtoes <- lapply(processnames[processnames != molecular_deposition_proces],
                           self$FromDataAndTo)
        otherairsinks <- fromtoes |>
          purrr::compact() |>
          purrr::keep(~ "from.SubCompart" %in% names(.x)) |>
          dplyr::bind_rows() |>
          dplyr::filter(from.SubCompart == "air") |>
          dplyr::pull(process) |>
          unique()

        # fromtoesflow <- lapply(flowNames,
        #                    self$FromDataAndTo)
        # names(fromtoesflow) <- flowNames
        # otherairflows <- fromtoesflow |>
        #   purrr::imap(~ dplyr::mutate(.x, flowname = as.character(.y))) |>
        #   purrr::compact() |>
        #   purrr::keep(~ "SubCompart" %in% names(.x)) |>
        #   dplyr::bind_rows() |>
        #   dplyr::filter(SubCompart == "air") |>
        #   dplyr::pull(flowname) |>
        #   unique()

        # create the airsinks reactive
        node_name <- "other_MolAir_sinks"
        
        other_MolAir_sinks <- make_combiner_fun(otherairsinks)
        
        # assign the function into the environment where it will be found
        assign(node_name, other_MolAir_sinks, envir = .GlobalEnv)  
        
        private$myReactiveDAG$add_node_reactive(
          node_name      = node_name,
          node_fun_name  = node_name
        )
      }
      
      function_names <- c(
        processnames, flowNames
      )
      
      # Initial needvars (params of the initial functions)
      needvars <- unlist(
          lapply(function_names, function(fn) private$get_clean_args(fn, ns)),
          use.names = FALSE
        ) |>
        unique()
      needvars <- needvars[!needvars %in% private$myReactiveDAG$knowns]
      missing_vars <- list()
      
      # complete SBvars
      repeat {
        needvars <- needvars[!needvars %in% c("parent", "...", "TRUE")] # parent resolved by ReactiveDAG, ... by R, TRUE??
        missing_vars <- c(missing_vars, needvars[!(needvars %in% my_ls | needvars %in% private$myReactiveDAG$knowns)])
        SBoofunction2do <- needvars[needvars %in% my_ls | needvars %in% private$myReactiveDAG$knowns]
        
        if (length(SBoofunction2do) == 0) {
          break
        }
        
        for (fun2do in SBoofunction2do) {
          private$myReactiveDAG$add_node_reactive(node_name = fun2do)
        }
        
        needvars <- unique(unlist(
          lapply(SBoofunction2do, function(fn) private$get_clean_args(fn, ns)$name),
          use.names = FALSE
          ))|>
            unique()
        needvars <- needvars[!needvars %in% private$myReactiveDAG$knowns]
        
      }
      
      # ---- determine order / complete / incomplete functions ----
      complete_nodes_piter <- list()
      known <- private$myReactiveDAG$dataNames
      
      testorder <- 0  
      repeat{
        testorder <- testorder + 1
        
        incomplete_nodes <- setdiff(private$myReactiveDAG$knowns, known)
        # which are now complete with known?
        complete_nodes <- character(0)
        for (fun in incomplete_nodes) {
          args_df <- private$get_clean_args(fun, ns)
          # optional_args <- args_df$name[args_df$has_default]
          # function is complete if all args are already known
          if (nrow(args_df) == 0 || all(args_df$name %in% known)) {
            complete_nodes <- c(complete_nodes, fun)
          }
        }
        if (length(complete_nodes) == 0){
          break
        }

        complete_nodes_piter[[testorder]] <- complete_nodes
        known <- c(known, complete_nodes)
      }
      
      # add flows and processes to DAG; collect air-processes
      for (pf in function_names){
        private$myReactiveDAG$add_node_reactive(node_name = pf)
        # # Molecular deposition depends on these
        # fromtos <- self$FromDataAndTo(pf) |>
        #   mutate(
        #     from_sub = coalesce(from.SubCompart, SubCompart)
        #   ) |>
        #   filter(Species == "Unbound", from_sub == "air")
        # if (nrow(fromtos) > 0 ){
        #   
        # }
      }

      # to debug vars, ...
      return (list(complete_nodes_piter = complete_nodes_piter, missing_vars = unique(missing_vars)))
    },  
    
    #' @description set a constant in the internal data, to enable use by SB variable etc.
    #' @param ... named value
    SetConst = function(...){
      #if ... is a list with a single data.frame, the const has dimensions
      Params <- list(...)
      #there can only be one at a time
      if (length(Params) != 1){
        stop("SetConst can deal with one Const at a time")
      }
      if (isa(Params, "list") && isa(Params[[1]], "data.frame")) {
        private$UpdateDL(Params[[1]])
      } else {
        private$UpdateDL(...)
      }
    },
    
    #' @description if UpdateRows is a list of variables (each containing a fetchdata() result): see private$MutateVar,
    #' if UpdateRows is (a csv filename of) a single dataframe it will be converted before the call
    mutateVars = function(UpdateRows) {
      #browser()
      if ("character" %in% class(UpdateRows)) {
        stopifnot(endsWith(UpdateRows, ".csv") && file.exists(UpdateRows))
        UpdateRows <- read.csv(UpdateRows) #now it's a data.frame
      }
      if ("data.frame" %in% class(UpdateRows)) {
        # you cannot test enough..
        if (!all(c("varName", "Waarde") %in% names(UpdateRows))) {
          stop("UpdateRows should contain columns varName and Waarde; dimensions when needed")
        }
        # fetch all variables, empty the dataframes and make room for new rows
        uniqvNames <- unique(UpdateRows$varName)
        baseVars <- lapply(uniqvNames, private$FetchData)
        names(baseVars) <- uniqvNames
        Updated <- lapply(baseVars, function(x) {
          if ("data.frame" %in% class(x)) {
            x[F,]
          } else {
            unlist(unname(x))
          }
        }) 
        names(Updated) <- names(baseVars)
        
        #put rows in Updated
        for (vari in 1:nrow(UpdateRows)){ #loop vnamesDistSD rows 
          
          vname <- UpdateRows$varName[vari]
          if ("data.frame" %in% class(baseVars[[vname]])) {
            #apply to new row of data to Mutate
            Dees <- The3D[The3D %in% names(UpdateRows)]
            NeedDees <- Dees[!is.na(UpdateRows[vari,Dees])]
            newRow <- UpdateRows[vari, NeedDees, drop = F]
            newRow[,vname] <- UpdateRows$Waarde[vari]
            Updated[[vname]] <- rbind(Updated[[vname]], newRow)
          } else { #just a number
            Updated[[vname]] <- UpdateRows$Waarde[vari] 
          }
        }  
        UpdateRows <- Updated
      } # now its a list of data.frames / names list value, like fetchdata() results
      
      for (i in 1:length(UpdateRows)) {
        if ("data.frame" %in% class(UpdateRows[[i]])) {
          private$MutateVar(UpdateRows[[i]])
        } else {
          UR <- UpdateRows[i] #mind the single []; leaving atomic values as list
          if (is.list(UpdateRows[[1]])) {
            UR[[1]] <- UR[[1]][[1]]
          }
          private$MutateVar(UR) 
        }
      }
    },
    
    #' @description Replaces a complete table in the internal data system. 
    #' Use with care.
    #' @param UpdateDF the new table, replcing the previous
    #' @param keys dims that determines the table to be updated, see whichDataTable
    #' @param TableName alternative to keys: provide the tablename. Use with care!
    UpdateData = function(UpdateDF, keys, TableName = NULL){
      stop("R6Core$UpdateData needs adapting to ReactiveDAG")
    },
    
    #' @description derive transfer processes from datalayer
    #' @param processName the one you are looking for, similar to  
    FromDataAndTo = function(processName){
      # Get all form-to state combination, given
      # [3D]Processes sheets and
      # Processes columns in [3D]Sheet
      # And restrict to existing states
      # k_Advection is a aggregation of all flows
      if (processName == "k_Advection"){
        warning("k_Advection is a aggregation of all flows")
        return(NULL)
      }
      hit3D <- vapply(
        The3D,
        function(d) any(self$RawTables$ProcessFromTo[[d]]$process == processName),
        logical(1)
      )
      if (sum(hit3D) == 0) { # it may be a flow!

        flows <- self$RawTables$FlowIO |> dplyr::filter(FlowName == processName)
        if (nrow(flows) == 0){
          stop(paste("no from/to found for", processName))
        }
        one3D <- unique(flows$Dimension)
        stopifnot(length(one3D) == 1)
        forwhichCol <- ifelse(one3D == "Scale", "SubCompart", "Scale") #the "other"
        dim1 <- flows |> dplyr::rename(
          !!paste("from", one3D, sep = ".") := from,
          !!paste("to",   one3D, sep = ".") := to,
          !!forwhichCol := forWhich
        ) |> 
          dplyr::select(-FlowName, -Dimension) |>
          expand.grid.df(self$states$asDataFrame |> 
                           dplyr::select(Species) |>
                           dplyr::distinct()
          )
        return(dim1)
      }
      # it is a process
      stopifnot(sum(hit3D) == 1) #there can only be one
      one3D <- The3D[hit3D]
      
      from.one3D <- paste("from", one3D, sep = ".")
      to.one3D <- paste("to", one3D, sep = ".")
      
      dim1 <- self$RawTables$ProcessFromTo[[one3D]] |>
        dplyr::filter(process == processName) |>
        dplyr::rename(
          !!from.one3D := from,
          !!to.one3D := to
        )
      
      ########## local helper 
      # dimension: "SubCompart" or "Species"
      # process: e.g. "k_DryDeposition"
      # values: vector of SubCompart or Species values to check
      is_allowed <- function(dimension, process, values) {
        exc <- self$RawTables$k_exceptions[[dimension]]
        
        # If process column not present in exception table: everything allowed
        if (!process %in% colnames(exc)) {
          return(rep(TRUE, length(values)))
        }
        
        # Match values against exception table and read the TRUE/FALSE flag
        idx <- match(values, exc[[dimension]])
        # TRUE/FALSE for those that exist in table; NA treated as TRUE (no rule)
        allowed <- exc[[process]][idx]
        allowed[is.na(allowed)] <- TRUE
        
        allowed
      }

      #expand to other dims
      Other <- The3D[!hit3D]
      
      for (oDim in Other){ 
        fulldim <- self$RawTables$DimsTables[[oDim]]
        allowed <- is_allowed(oDim, processName, fulldim[[oDim]])
        dim1 <- expand.grid.df(dim1, fulldim[allowed, oDim])
        names(dim1)[names(dim1) == "y"] <- oDim
      }
      
      # test if both from and to states present in self$states
      dimfrom <- dim1 |>
        dplyr::mutate(.row_id = dplyr::row_number()) |> #keep rownumbers to keep from dim1
        dplyr::select(Other, !!from.one3D, .row_id) |>
        dplyr::rename(!!one3D := !!from.one3D) |>
        self$states$clipStates() 
      dimto <- dim1 |>
        dplyr::mutate(.row_id = dplyr::row_number()) |> #keep rownumbers to keep from dim1
        dplyr::select(Other, !!to.one3D, .row_id) |>
        dplyr::rename(!!one3D := !!to.one3D) |>
        self$states$clipStates() 
      # delete those rows for which either from-to state do not exist
      if (nrow(dimfrom) != nrow(dim1) | nrow(dimto) != nrow(dim1)) {
        existrows <- union(
          dimfrom |> dplyr::pull(.row_id),
          dimto |> dplyr::pull(.row_id)
        )
        dim1 <- dim1 |>
          dplyr::slice(existrows)
      }
      return(dim1)
    },
    
    #' @description states and all data-layer to an rds-file
    save_world = function(filename){
      whole_world <- list(
        substance = private$Substance,
        states = private$States$asDataFrame,
        DAGdata = private$myreactiveDAG$sources) # lapply fetchdata?
        # private$myreactiveDAG$params ? not (yet) in use
        # not used anymore filterstates = private$filterstates)
      saveRDS(whole_world, file = filename)
    },
    
    load_world = function(filename){
      whole_world <- readRDS(filename)
      private$Substance <- whole_world$substance
      private$States <- SBstates$new(whole_world$states)
      private$States$myCore <- self
      private$myreactiveDAG$sources <- whole_world$DAGdata 
      # not used anymore private$filterstates <- whole_world$filterstates
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
    Flows = function(value){
      if (missing(value)) {
        private$flows
      } else {
        stop("`$Flows` are set at initialise", call. = FALSE)
      }
    },
    Processes = function(value){
      if (missing(value)) {
        private$processes
      } else {
        stop("`$Processes` are set at initialise", call. = FALSE)
      }
    },
    substancemessage = function(value){
      if (missing(value)) {
        private$Substancemessage
      } else {
        private$Substancemessage <- value
      }
    },
    #' @field kaas getter for r.o. property (all k's)
    kaas = function(value) {
      if (missing(value)) {
        if (is.null(private$SBkaas)) return(NULL)
        #else
        private$SBkaas
      } else {
        warning("`$kaas` should be set by UpdateKaas(), format is strict!", call. = FALSE)
        private$SBkaas <- value
      }
    },
    
    SuBmode = function(value){
      if (missing(value)){
        private$suBmode
      } else {
        if (value %in% c("Molecular", "Particle")){
          private$SBmode <- value
        } else {
          # should be a list or path to jsonfile containing
          private$suBmode <- if (is.list(value)) {
            value
          } else {
            stopifnot(file.exists(file.path(private$myDataLocation, value)))
            jsonlite::fromJSON(file.path(private$myDataLocation, value), simplifyDataFrame = FALSE)
          }
          # extract needed SBmode
          private$SBmode <- {
            SubstanceList <- private$suBmode[["Substance"]]
            SubstanceList$SBmode
          }
          
        }
      }
    },
    substance = function(value){
      if (missing(value)){
        private$Substance
      } else {
        # update Substance and other data from 
        # substance block in private$suBmode to ReactiveDAG
        theList <- private$suBmode[[value]]
        browser() # are there multiple substances??
        private$set_sources(theList)
        private$Substance <- value
      }
    },
    RawTables = function(value){
      if (missing(value)){
        private$rawTables
      } else {
        # update Substance in ReactiveDAG
        private$rawTables <- value
      }
    },
    RowIdentifyers = function(value){
      if (missing(value)){
        private$rowIdentifyers
      } else {
        # update Substance in ReactiveDAG
        private$rowIdentifyers <- value
      }
    },

    #' @description Export everything in fetchdata
    exportMetadata = function() {
      #browser()
      fetchdatanames <- unique(self$fetchData())

      exclusions <- c("Dimension", "FlowName", "forWhich", "Unit", "VarName", "X", "DefaultFRACarea", "DefaultPH", "Kow_default", "MOLMASSAIR", "outdated", "outdated.1", "outdated.2", "Pvap25_default")
      
      # Initialize an empty tibble to store the results
      result_tibble <- tibble(
        varname = character(),  # Column for the variable name
        data = list()           # Column for nested tibbles
      )
      
      # Iterate through the fetched data names
      for(i in fetchdatanames) {
        if (!i %in% exclusions) {
          
          fd <- data.frame(self$fetchDataUnits(i))
          
          # Convert fd to a tibble if it's not already one
          fd_tibble <- as_tibble(fd)
          
          # Add a new row to the result_tibble
          result_tibble <- bind_rows(result_tibble, tibble(
            varname = i,       # Column for the variable name
            data = list(fd_tibble)  # Nest the tibble 
          ))
        }
      }
      
      result_tibble <- result_tibble |>
        distinct()
      
      # Get the kaas and format it the same as the fetchdata items
      kaas <- self$kaas
      kaas_tibble <- as_tibble(kaas)
      kaas_tibble <- tibble(varname = "kaas",
                            data = list(kaas))
      
      # Bind the fetchdata tibble and the kaas tibble together
      result_tibble <- rbind(result_tibble, kaas_tibble)
    },
    
    filterStates = function(value){
      stop("use methods of SBstates class")
    }
  ),
  
  private = list(
    myDataLocation = NULL,
    SBmode = NULL,
    suBmode = NULL,
    debugR = NULL,
    rawTables = NULL,
    rowIdentifyers = NULL,
    myReactiveDAG = NULL,
    
    States = NULL,
    Substance = NULL,
    Substancemessage = FALSE,
    solvername = NULL, 
    SBkaas = NULL,
    ModuleList =  NULL,
    nodeList = NULL,
    solver = NULL,
    l_postPoneList = NULL,
    concentration = NULL,
    flows = NULL,
    processes = NULL,
    filterstates = list(),
    substanceproperties = list(),

    
    getStates = function(Scales, Species, SubComparts,
                         Rules, ModeRules){

      # create States -- MyCore Needed ??
      states = SBstates$new(Scales, SubComparts, Species)
      
      states$filterStates(Rules)
      states$filterStatesByMode(private$SBmode, ModeRules)
      
      return(states)
    },
    
    # Get block from json structure
    Block2DAG = function(theList){
      
      purrr::imap(theList, function(aVar, nm) {
        # Lege waardes bij pKa bijv werden door bind_rows een lege tibble, extra safeguard gebouwd
        # if (nm =='pKa') {
        #   browser()
        # }
        # 
        # if (!is.atomic(aVar)) {
        #   aVar <- aVar |>
        #     dplyr::bind_rows(.id = "id")
        # }
        if (is.null(aVar)) {
          aVar <- NA
        } else if (is.list(aVar) && length(aVar) == 0) {
          aVar <- NA
        } else if (!is.atomic(aVar)) {
          aVar <- aVar |>
            dplyr::bind_rows(.id = "id")
        }
        
        # Apply SI conversion if exists
        if (!is.atomic(aVar) && !is.null(self$RawTables$Units)) {
          SIexpression <- self$RawTables$Units$ToSI[
            self$RawTables$Units$VarName == nm
          ]
          if (length(SIexpression) == 1 && !is.na(SIexpression) && SIexpression != "") {
            # Apply conversion to the data frame column
            if (nm %in% names(aVar)) {
              aVar[[nm]] <- eval(parse(text = gsub(nm, "aVar[[nm]]", SIexpression)))
            }
          }
        } else if (is.atomic(aVar) && !is.null(self$RawTables$Units)) {
          # For scalar values
          SIexpression <- self$RawTables$Units$ToSI[
            self$RawTables$Units$VarName == nm
          ]
          if (length(SIexpression) == 1 && !is.na(SIexpression) && SIexpression != "") {
            assign(nm, aVar)
            aVar <- eval(parse(text = SIexpression))
          }
        }
        
        private$myReactiveDAG$set_source(nm, aVar)
      })
    },
    
    # helper: for a single value (or vector), try numeric, then logical, else keep as is
    coerce_scalar = function(x) {
      # assume length-1 here in vector mode; still works for length > 1
      num <- suppressWarnings(as.numeric(x))
      if (!any(is.na(num))) {
        return(num)
      }
      
      logi <- suppressWarnings(as.logical(x))
      if (!any(is.na(logi))) {
        return(logi)
      }
      
      x
    },
    
    # use data from Defs, possibly towider, .. 
    # and put each var into a reactive value, 
    # possibly inherit Matrix, Component
    data2DAG = function(DefsTable, Table_type = "vector") {
      
      Doexpression <- function (varname, x, SIexpression){ #execute expression to convert to SI
        #NB varname is local here
        assign(varname, x)
        eval(parse(text = SIexpression))
      }
      
      dims <- self$RawTables$rowIdentifyers
      dim2use <- dims[dims %in% names(DefsTable)]
      
      if ("Substance" %in% dim2use) {
        if (!self$substancemessage) {
          message("Substance properties now in separated json-file / list()")
          self$substancemessage <- TRUE
        }
        return()
      }
      
      switch(
        Table_type,
        
        # default: each non-dimension column becomes a scalar/vector source
        "vector" = {
          for (nam in names(DefsTable)){
            aValue <- private$coerce_scalar(unlist(DefsTable[nam]))
            
            #Convert to SI
            if (!is.null(self$RawTables$Units)) {
              SIexpression <- self$RawTables$Units$ToSI[
                self$RawTables$Units$VarName == nam
              ]
              if (length(SIexpression) == 1 && !is.na(SIexpression) && SIexpression != "") {
                aValue <- Doexpression(nam, aValue, SIexpression)
              }
            }
            
            private$myReactiveDAG$add_source(nam, aValue)
          }
        },
        
        # variables in columns
        "Vars" = {
          
          if (length(dim2use) == 3 &&
              all(dim2use %in% c("SubCompart", "Matrix", "Compartment"))) {
            dim2use <- "SubCompart"
          }
          
          cols_to_use <- names(DefsTable)[!names(DefsTable) %in% dim2use]
          # exclude those with a process name (indicating from and to, not reactive!)
          cols_to_use <- cols_to_use[!startsWith(cols_to_use, "k_")]
          
          tibble_list <- purrr::map(
            cols_to_use,
            ~ DefsTable |> dplyr::select(dplyr::all_of(c(dim2use, .x)))
          )
          names(tibble_list) <- cols_to_use
          
          tibble_list <- tibble_list |>
            purrr::imap(function(.x, .y) {
              
              # Convert to SI units
              if (!is.null(self$RawTables$Units)) {
                SIexpression <- self$RawTables$Units$ToSI[
                  self$RawTables$Units$VarName == .y
                ]
                if (length(SIexpression) == 1 && !is.na(SIexpression) && SIexpression != "") {
                  .x[,.y] <- Doexpression(.y, .x[,.y], SIexpression)
                }
              }
              .x
            })
          tibble_list |>
            purrr::iwalk(~ private$myReactiveDAG$add_source(.y, .x))
        },
        
        # from long format, variables in VarName, Waarde
        "ToWiden" = {
          needcols <- c(dim2use, "VarName", "Waarde")
          
          tibble_list <- DefsTable |>
            dplyr::filter(!is.na(Waarde)) |>
            dplyr::select(dplyr::all_of(needcols)) |>
            split(~ VarName)
          
          for (nm in names(tibble_list)) {
            
            # Convert to SI-units
            if (!is.null(self$RawTables$Units)) {
              SIexpression <- self$RawTables$Units$ToSI[
                self$RawTables$Units$VarName == nm
              ]
              if (length(SIexpression) == 1 && !is.na(SIexpression) && SIexpression != "") {
                tibble_list[[nm]] <- tibble_list[[nm]] |> 
                  dplyr::mutate(Waarde = Doexpression(nm, Waarde, SIexpression))
              }
            }
            
            exclVarNam <- tibble_list[[nm]] |>
              dplyr::select(-VarName) |>
              dplyr::rename(!!nm := Waarde)
             
            private$myReactiveDAG$add_source(nm, exclVarNam)
          }
        }
      )
    },
    
    get_clean_args = function(fun_name, ns = getNamespace("sboo")) {
      f    <- get(fun_name, envir = ns)
      fmls <- formals(f)
      
      arg_names  <- names(fmls)
      args_clean <- gsub("^(to|from|all)\\.", "", arg_names)
      
      has_default <- vapply(fmls, function(x) !is.null(x), logical(1))
      
      data.frame(
        name        = args_clean,
        has_default = has_default,
        stringsAsFactors = FALSE
      )
    },
    
    # list of a vector and dataframes to be moved to private$myReactiveDAG
    set_sources = function(aList){
      browser()
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
    
    PrepareKaaslist = function(kaaslist){
      # find columnname starting with k_
      # rename to "k" and filter values > 0 for that column
      # find changing dim for this process, extend with from, to for other dim
      kaaslist <- kaaslist[!sapply(kaaslist, anyNA)]
      kaaslist <- lapply(kaaslist, function (process){
        
        process.names <- names(process)
        k.name <- process.names[startsWith(process.names,"k_")]
        otherdimNames <- process.names[
          !(process.names %in% c("process", k.name) |
              startsWith(process.names, "to.") |
              startsWith(process.names, "from.")
          )
        ]
        for (o1 in otherdimNames){
          to.o1 <- paste0("to.", o1)
          from.o1 <- paste0("from.", o1)
          process <- process |>
            dplyr::mutate(!!from.o1 := .data[[o1]]) |>
            dplyr::rename(!!to.o1 := !!o1)
        }
        # return
        process |>
          dplyr::rename(k = !!k.name)
      })
      VolleKaas <- do.call(rbind, kaaslist)
      
      if(length(VolleKaas) > 0) { #length is 0 if no processes 
        return(VolleKaas)
      } else {return(NULL)}
      
    },
    

    #' @description update rows of a variable in the "datalayer" (SB4Ndata-dataframes)
    #' @param UpdateRows the data.frame with the variable, dimensions should match and  
    #' 1 other column containing values, named as the var name
    #' @return side-effect; new "fetchdata()" is returned
    MutateVar = function(UpdateRows) {
      browser()
      if ("data.frame" %in% class(UpdateRows)) {
        dims <- The3D[The3D %in% names(UpdateRows)]
        varName <- names(UpdateRows)[!names(UpdateRows) %in% dims]
        if (length(varName) != 1){
          stop(paste("illegal call to MutateVar; multiple variable names found:", varName))
        }
        nowVar <- private$FetchData(varName)
        if (anyNA(nowVar)) {
          stop(paste("non-existing variable:", varName))
        }
        if (varName %in% names(private$moduleList)) {
          stop(paste("illegal call to MutateVar; variable is a calculated value:", varName))
        }
        #merge the "newValue" to the nowVar
        names(UpdateRows)[names(UpdateRows) == varName] <- "newValue"
        wannebeVar <- merge(nowVar, UpdateRows, all.x = T) |>
          distinct()
        if (nrow(nowVar) < nrow(wannebeVar)) {
          #browser() #start debugging?
          stop("illegal UpdateRows in MutateVar")
        }
        #update wannebe, remove newValue
        wannebeVar[!is.na(wannebeVar$newValue), varName] <- wannebeVar$newValue[!is.na(wannebeVar$newValue)]
        wannebeVar$newValue <- NULL
        private$myReactiveDAG$add_source()
        
      } else {# just a number
        private$myReactiveDAG$add_source()
      }
    },

    #' @description calculate and update a variable in the "datalayer" (SB4Ndata-dataframes)
    #' @param VarFunName the name of the variable to be updated AND the name of the function
    #' @param DIMRestrict optional list of restrictions; each element is c(DimName, Comparator, Value)
    #' uh, not in use at this point in time
    #' return side-effect; vector?
    UpdateDL = function(VarFunName = NULL, DIMRestrict = NULL, ...) {
      stop("use the methods related to reactiveDAG")
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
