#' @description
#' A short description...
#' 
#' @import R6
#' @export
EmissionModule <-
  R6::R6Class(
    "EmissionModule",
    public = list(
      initialize = function(TheSolver, input, ...) {#Solver is reference to States
        private$MySolver <- TheSolver
        MoreParams <- list(...)
        #switch between option 1) read from csv 2) read from excel 
                # 3) list of functions or 4) a data.frame 
        if (length(MoreParams) < 1){
          stop("expected filename or list of functions for emissions")
        }
        if ("character" %in% class(input)) {
          switch (tools::file_ext(input,
            "csv" = private$readfromcsv(input, ...),
            "xlsx" = private$readFromExcel(input, ...)
          )) 
        } else {# [[1]] a list of functions, named by States
          if ("list" %in% class(input && "function" %in% input[[1]])) {
            private$setEmissionFunctions(input)
          } else {
            private$setEmissionDataFrame(input)
          }
        }
        if ("unitFactor" %in% names(MoreParams)){
          private$UnitFactor <- MoreParams[["unitFactor"]]
        }
      },
      
      emissions = function(scenario = NULL){
        if (is.NULL(private$emissions)) return (NULL)
        if (is.null(scenario)) 
          return (names(private$emissions))[!names(private$emissions) %in% c(
                                           "Scenario","VarName","Unit")]
      }
      
    ),

    private = list(
      MySolver = NULL,
      Emissions = NULL,
      EmissionSource = NULL,
      UnitFactor = 1,
      Scenarios = NULL,
      Times = NULL,
      readFromClassicExcel = function(fn) {
        tryCatch(df <- openxlsx::read.xlsx(fn, sheet = "scenarios", startRow = 3),
                 error = function(e) NULL)
        if (is.null(df)) stop ("Not a proper scenarios sheet in the xlsx")
        #clean a lot
        names(df)[c(3,4)] <- c("VarName", "Unit") #names(df) == "X3"
        df$current.settings <- NULL
        #remove all Xnames
        Xnames <- startsWith(names(df), "X")
        df <- df[startsWith(df$VarName, "E.") & !is.na(df$VarName),!Xnames]
        StateAbbr <- substr(df$VarName, start = 3, stop = 9)
        df$i <- private$MyCore$findState(StateAbbr)
        df
      },
      
      readfromexcel = function(fn, ...){
        sheetNames <- openxlsx::getSheetNames(fn)
        #ignore the Sheetx names others are assumes to be scenario name
        sheetNames <- sheetNames[!grepl("Sheet[23]", sheetNames)]
        if (length(sheetNames) == 1){
          private$EmissionSource <- openxlsx::read.xlsx(fn, sheet = sheetNames)
        } else {
          dfs <- lapply(sheetNames, function(sheet){
            openxlsx::read.xlsx(fn, sheet = sheet)
          })
          for (i in 1:length(sheetNames)){
            df <- dfs[[i]]
            df$scenario <- sheetNames
          }
          private$EmissionSource <- do.call(rbind, dfs)
        }
      },
      
      readfromcsv = function(fn, ...) {
        df <- read.csv(fn)
        setEmissionDataFrame(df)
      },
      
      setEmissionDataFrame = function(emis_df) {
        if ("data.frame" %in% class(emis_df) && all(c("Abbr","Emis") %in% names(emis_df))) {
          #we need states - via solver from the core
          states <- private$MySolver$solveStates
          # there can be multiple times in optional Timed column. 
          #if so, emis_df becomes a list of vectors
          if ("Timed" %in% names(emis_df)) {
            Times <- sort(unique(emis_df$Timed))
            vEmis <- replicate(states$nStates, data.frame(Timed = Times[1], emis = 0.0), simplify = F)
            names(vEmis) <- states$asDataFrame$Abbr
            #update the first time if present, 
            #then append all emis_df in the right state list
            #make sure the last Times is present for all states; set to 0.0 if missing
            timedRows <- emis_df[emis_df$Timed == Times[1],]
            for (irow in 1:nrow(timedRows)){
              vEmis[[timedRows$Abbr[irow]]]$emis <- timedRows$Emis[irow] # 1000000 / Molweight / (3600*24*365) #t/an -> mol/s
            }
            for (atime in Times[2:(length(Times)-1)]) {
              timedRows <- emis_df[emis_df$Timed == atime,]
              for (irow in 1:nrow(timedRows)){
                newRow <- data.frame(Timed = timedRows$Timed[irow], 
                                     emis = timedRows$Emis[irow])  # 1000000 / Molweight / (3600*24*365)) #t/an -> mol/s)
                vEmis[[timedRows$Abbr[irow]]] <- rbind(vEmis[[timedRows$Abbr[irow]]], newRow)
              }
            }
            lastTime <- tail(Times,1)
            timedRows <- emis_df[emis_df$Timed == lastTime,]
            for (aState in names(vEmis)) {
              posEmis <- timedRows$Emis[timedRows$Abbr == aState]  # 1000000 / Molweight / (3600*24*365) #t/an -> mol/s
              emisMust <- ifelse(length(posEmis) > 0, posEmis, 0.0)
              newRow <- data.frame(Timed = lastTime, emis = emisMust)
              vEmis[[aState]] <- rbind(vEmis[[aState]], newRow)
            }
            private$EmissionSource <- vEmis
          } else { 
            #steady state
            vEmis <- rep(0.0, length.out = states$nStates)
            names(vEmis) <- states$AsDataFrame$Abbr
            vEmis[match(emissions$Abbr, states$asDataFrame$Abbr)] <- emissions$Emis
            # from kg/yr to Mol/s
            private$EmissionSource <- vEmis# * 1000000 / Molweight / (3600*24*365) #t/an -> mol/s
            names(private$EmissionSource) <- states$asDataFrame$Abbr
            
          }
        }
      }
    )
  )