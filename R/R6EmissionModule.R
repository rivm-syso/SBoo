#' @title Emission module
#' @description to set/get emissions, possibly make into functions
#' @import R6
#' @export
EmissionModule <-
  R6::R6Class(
    "EmissionModule",
    public = list(
      initialize = function(mySolver, emis, uncertainFun = list()) {
        #browser()

        private$mySolver = mySolver
        private$uncertainFun <- uncertainFun
        # determine type (vector / dataframe, with timed, dynamic(as dataframe to convert or list of function))
        if ("data.frame" %in% class(emis)) {
          private$setEmissionDataFrame(emis)
                    
          #character ? -> read from file as vector
        } else if (is.character(emis)){
          switch (tools::file_ext(emis,
               "csv" = private$readfromcsv(emis),
               "xlsx" = private$readFromExcel(emis)
          )) 
          private$emission_tp <- "Vector"
          
          # list of functions?
        } else if ("list" %in% class(emis)) {
          for(i in emis){
            if(!is.function(i)){
              stop("Not all elements in the list are functions")
            }
          }
          private$emission_tp <- "DynamicFun"
          private$Emissions <- emis          
        }
        
      },
      
      
      emissions = function(scenario = NULL){ #scenario also used for "RUN"
        if (is.null(private$Emissions)) return (NULL)
        if (!is.null(scenario)) 
          return (names(private$Emissions))[!names(private$Emissions) %in% c(
            "Scenario","VarName","Unit")]

        if (private$emission_tp != "Vector") {
          stop("emissions cannot be casted to a named vector")
        }
        #normal use? expand to all state in mySolver and return as named vector
        if (is.na(RUNs)) {
          emis <- dplyr::left_join(private$mySolver$solveStates, private$Emissions)
          emis$Emis[is.na(emis$Emis)] <- 0
          return(emis)
        } else { #apply uncertainty RUNs
          # RUNs should contain lhs (like) 0-1 numbers per relevant state
          
        }
      },
      
      # return approx function 
      emissionFunctions = function(states) {
        
          if (private$emission_tp == "DynamicFun"){
            if (!all(names(private$Emissions) %in% states$asDataFrame$Abbr)) {
              notfound <- as.list(names(private$Emissions)[
                !names(private$Emissions) %in% states$asDataFrame$Abbr])
              stop(do.call(paste,c(list("not all states in SB engine (the matrix)"), notfound)))
            }
            if (! is.na(private$uncertainFun)){
              stop("not possible to combine uncertain emissions with dynamic emissions, yet")
            }
            return(private$Emissions)
          }
        
          if (private$emission_tp == "Dynamic_df") {
            #return the df as list of functions
            if(!(all(c("Abbr","Emis", "Timed") %in% names(private$Emissions)))){
              stop("Expected 'Abbr', 'Emis' and 'Timed' columns in dataframe")
            }
            if (! is.na(private$uncertainFun)){
              stop("not possible to combine uncertain emissions with dynamic emissions, yet")
            }
            
            if(!all(as.character(states$asDataFrame$Abbr) %in% as.character(states$asDataFrame$Abbr))){
              stop("Abbreviations are not compatible with states")
            }
            #make 'm
            return(private$makeApprox(private$Emissions))
          }
          # else
          stop("no dynamic emissions") #or make them a level line ???
      }
      
    ),

    private = list(
      mySolver = NULL,
      emission_tp = NULL,
      Emissions = NULL, #vector / dataframe or list of functions, attributes as input at init
      uncertainFun = NA,
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
        setEmissionDataFrame(df)
      },
      
      readfromexcel = function(fn){
        sheetNames <- openxlsx::getSheetNames(fn)
        #ignore the Sheetx names others are assumes to be scenario name
        sheetNames <- sheetNames[!grepl("Sheet[23]", sheetNames)]
        if (length(sheetNames) == 1){
          private$setEmissionDataFrame(openxlsx::read.xlsx(fn, sheet = sheetNames))
        } else {
          dfs <- lapply(sheetNames, function(sheet){
            openxlsx::read.xlsx(fn, sheet = sheet)
          })
          for (i in 1:length(sheetNames)){
            df <- dfs[[i]]
            df$scenario <- sheetNames
          }
          private$setEmissionDataFrame(do.call(rbind, dfs))
        }
      },
      
      readfromcsv = function(fn) {
        df <- read.csv(fn)
        setEmissionDataFrame(df)
      },

      setEmissionDataFrame = function(emis) {
        if (all(c("Abbr","Emis") %in% names(emis))) {
          if ("Timed" %in% names(emis)){
            private$emission_tp <- "Dynamic_df"
          }
          else {
            private$emission_tp <- "Vector"
          }
        } else stop ("at least columns with names Abbr,Emis")
        private$Emissions <- emis
      },
      
      # Create function to make approx functions from data (input is a df with the columns Abbr, Timed and Emis)
      makeApprox = function(vEmissions, states){
        is.df.with(vEmissions, "EmissionModule$makeApprox", c("Timed", "Emis", "Abbr"))
        
        vEmis <- 
          vEmissions |> 
          group_by(Abbr) |> 
          summarise(n=n(),
                    EmisFun = list(
                      approxfun(
                        data.frame(Timed = c(0,Timed), 
                                   Emis=c(0,Emis)),
                        rule = 1:2)
                    )
          )
        funlist <- vEmis$EmisFun
        names(funlist) <- vEmis$Abbr
        return(funlist)
      } 
      
    )
  )