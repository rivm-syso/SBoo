# ReadRawData, just the sheets
ReadRawData <- R6::R6Class("SBreadData",
  public = list(
    initialize = function(SBooDataLocation){
      
      stopifnot(dir.exists(SBooDataLocation))
      # site effect create private ToWidenList, VarsList, OtherList
      private$readMasCsvs(SBooDataLocation)  
      
  }),
  
  active = list(
    rowIdentifyers = function(value){
      if (missing(value)) {
        private$RowIdentifyers
      } else {
        stop("`$rowIdentifyers` are R.O>", call. = FALSE)
      }
    },
    tableDims = function(value){
      if (missing(value)) {
        c(The3D, "Matrix", "Compartment")
      } else {
        stop("`$tableDims` are R.O>", call. = FALSE)
      }
    },
    DimsTables = function(value){
      if (missing(value)) {
        private$DimList
      } else {
        stop("`$DimsTables` are R.O>", call. = FALSE)
      }
    },
    VarsTables = function(value){
      if (missing(value)) {
        private$VarsList
      } else {
        stop("`$VarsTables` are R.O>", call. = FALSE)
      }
    },
    ToWidenTables = function(value){
      if (missing(value)) {
        private$ToWidenList
      } else {
        stop("`$ToWidenTables` are R.O>", call. = FALSE)
      }
    },
    SubstanceProperties = function(value){
      if (missing(value)) {
        private$substanceProperties
      } else {
        stop("`$SubstanceTables` are R.O>", call. = FALSE)
      }
    },
    CONST = function(value){
      if (missing(value)) {
        private$CONSTANTS
      } else {
        stop("`$SubstanceTables` are R.O>", call. = FALSE)
      }
    },
    ProcessFromTo = function(value){
      if (missing(value)) {
        private$processfromto
      } else {
        stop("`$processfromto` are R.O>", call. = FALSE)
      }
    },
    
    FlowIO = function(value){
      if (missing(value)) {
        private$flowIO
      } else {
        stop("`$FlowIO` are R.O>", call. = FALSE)
      }
    },
      
    Units = function(value){
      if (missing(value)) {
        private$units
      } else {
        stop("`$units` are R.O>", call. = FALSE)
      }
    },
    k_exceptions = function(value){
      if (missing(value)) {
        private$K_exceptions
      } else {
        stop("`$k_exceptions` are R.O>", call. = FALSE)
      }
    }
  ),
  
  private = list(
    ToWidenList = list(), 
    VarsList = list(),
    DimList = list(),
    K_exceptions = list(),
    SubstanceList = list(),
    processfromto = list(),
    CONSTANTS = NULL,
    units = NULL,
    flowIO = NULL,
    
    RowIdentifyers = c(
      The3D,
      paste("to", The3D, sep = "."),
      "Substance", "process", "Compartment", "Matrix", "from", "to", "QSAR.ChemClass"
    ),
    
    readMasCsvs = function(MlikeFile){

      OtherDimColNames = c("Default", "AbbrS", "AbbrC", "AbbrP",
                            paste(The3D, "Order", sep = ""))
      dimNames = paste(The3D, "Name", sep = "")
      
      ToWidenList <- list()
      VarsList <- list()
      DimList <- list()
      ProcessList <- list()
      k_exceptions <- list()
      substanceProperties <- list()

      for (Def in Defs) {
        tableName <- read.csv(
          paste(MlikeFile, "/", Def, ".csv", sep = ""))
        if (Def == "CONSTANTS") {
          private$CONSTANTS <- tableName$Waarde
          names(private$CONSTANTS) <- tableName$VarName
          next
        }
        if (Def == "FlowIO") {
          private$flowIO <- tableName
        }
        if (Def == "Units") {
          private$units <- tableName
          next
        }
        
        if ("VarName" %in% names(tableName)) {
          DimNames <- names(tableName)[names(tableName) %in% private$RowIdentifyers]
          ToWidenList[[length(ToWidenList)+1]] <- tableName[ ,c(DimNames, "VarName", "Waarde")]
          next
        } else {
          if (grepl("Processes", Def)) {
            procesDim <- sub("(.*)Processes.*", "\\1", Def)
            stopifnot(procesDim %in% The3D)
            ProcessList[[procesDim]] <- tableName
            next
          } else {
            wantColumns <- names(tableName)[!names(tableName) %in% OtherDimColNames] 
            VarsList[[length(VarsList)+1]] <- tableName[ ,wantColumns]
          }
        }
        Dims <- names(tableName)[names(tableName) %in% self$tableDims] 
        if (length(Dims) == 1 | (length(Dims) == 3 && all(Dims %in% c("SubCompart","Compartment","Matrix")))){
          # nasty exception; "Compartment", "Matrix" are properties here, for left_join()
          if (length(Dims) == 3) {
            TheDim = "SubCompart"
          } else {
            TheDim = Dims[1]
          }
          # separate from-to info from other
          k_columns <- startsWith(names(tableName), prefix = "k_")
          k_columns <- names(tableName)[k_columns]
          if (length(k_columns) > 0 ){
            k_exceptions[[TheDim]] <- tableName |> dplyr::select(all_of(c(TheDim, k_columns)))
            tableName <- tableName |> dplyr::select(-all_of(k_columns))
          }
          wantColumns <- c(Dims, names(tableName)[names(tableName) %in% OtherDimColNames]) 
          DimList[[TheDim]] = tableName[, wantColumns]
        }
      }
      
      private$ToWidenList <- ToWidenList 
      private$VarsList <- VarsList
      private$DimList <- DimList
      private$K_exceptions <- k_exceptions
      private$processfromto <- ProcessList
    }
  )
)