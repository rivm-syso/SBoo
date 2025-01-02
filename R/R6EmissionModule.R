#' @title Emission module
#' @description to set/get emissions, possibly make into functions
#' @import R6
#' @export
EmissionModule <-
  R6::R6Class(
    "EmissionModule",
    public = list(
      initialize = function(emis, solvedAbbr) {

        private$solvedAbbr <- solvedAbbr # if NULL it should be given with emis
        
        if (is.character(emis)){ # read from file as data frame
          emis <- switch (tools::file_ext(emis,
                                  "csv" = private$readfromcsv(emis),
                                  "xlsx" = private$readFromExcel(emis)
          )) 
        }
        
        if ("matrix" %in% class(emis)) {
          stopifnot(all(colnames(emis) %in% private$solvedAbbr))
          private$emission_tp <- "runs_row"
        } else {
          if ("list" %in% class(emis) & all(is.function(emis))) {
            private$emission_tp <- "DynamicFun"
          } else {
            if ("data.frame" %in% class(emis)){
              private$emission_tp <- private$setEmissionDataFrame(emis)
            } else stop ("unknown format of emissions")
          }
          
        }
        private$Emissions <- emis          
      },
      
      emissions = function(scenario = NULL){ #scenario also used for "RUN"
        
        if (is.null(private$Emissions)) return (NULL)

        if (!private$emission_tp %in% c("Vector", "runs_long", "runs_row")) {
          stop("emissions cannot be casted to a named vector")
        } 
        
        #normal use or scenarios / uncertainty runs
        emis <- rep(0, length(private$solvedAbbr))
        names(emis) <- private$solvedAbbr
          
        if (is.null(scenario)) {
          if (private$emission_tp != "Vector")
            stop("scenario/run is missing in emission data")
          else {
            emis[private$Emissions$Abbr] <- private$Emissions$Emis
            return(emis)
          }
        } 

        # RUNs should contain lhs (like) samples per relevant state rowwise or in long format
        if (private$emission_tp == "runs_row"){
          if (is.numeric(scenario)) {
            rowNum <- scenario
          } else {
            rowNum <- match(scenario, rownames(private$Emissions))
            stopifnot(is.na(rowNum))
          }
          lhssamples <- private$Emissions[rowNum,]
          emis[names(lhssamples)] <- lhssamples
        } else {
          lhssamples <- private$Emissions$Emis[private$emis$run == scenario,]
          emis[lhssamples$Abbr] <- lhssamples$Emis
        }
        
        return(emis)
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
      emission_tp = NULL,
      Emissions = NULL, #vector / dataframe or list of functions, attributes as input at init
      solvedAbbr = NULL, #vector Abbr of solveStates
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
          return(openxlsx::read.xlsx(fn, sheet = sheetNames))
        } else {
          dfs <- lapply(sheetNames, function(sheet){
            openxlsx::read.xlsx(fn, sheet = sheet)
          })
          for (i in 1:length(sheetNames)){
            df <- dfs[[i]]
            df$scenario <- sheetNames
          }
          return(do.call(rbind, dfs))
        }
      },
      
      readfromcsv = function(fn) {
        return (read.csv(fn))
      },

      # detemine format, return emission_tp
      setEmissionDataFrame = function(emis) {
        browser()
        if (all(c("Abbr", "Emis") %in% names(emis))) {
          if ("Timed" %in% names(emis)) {
            return("Dynamic_df")
          } else {
            stopifnot(all(emis$Abbr %in% private$solvedAbbr))
            if (any(c("run", "scenario") %in% names(emis))) {
              return ("runs_long")
            } else {
              return("Vector")
            }
          }
        } else {
          # it can be a data.frame with a run/scenario per row
          if (all(names(emis) %in% private$solvedAbbr)) {
            return("runs_row")
            
          } else stop ("at least columns with names Abbr,Emis or names equal to Abbr of states")
        }
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