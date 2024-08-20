#' @title Emission module
#' @description Function to get emissions from excel/csv files
#' @import R6
#' @export
EmissionModule <-
  R6::R6Class(
    "EmissionModule",
    public = list(
      initialize = function(input, ...) {#Solver is reference to States
        #browser()
        MoreParams <- list(...)
        #switch between option 1) read from csv 2) read from excel 
                # 3) list of functions or 4) a data.frame 
        
        emis <- MoreParams[[1]] # Get the emissions
        SF <- MoreParams[[2]] # Get the solver function 
        SB.K <- MoreParams[[3]] # Get the kaas 
        
        # First check if approxfuns are used in solver
        if("ApproxFun" %in% names(formals(SF))){
          private$setEmissionFunction(emis, SB.K)
        } 
        
        # else, check if the input is a data frame
        else if (class(emis) == "data.frame"){
          #ip <- Filter(is.data.frame, MoreParams)[[1]]
          private$setEmissionDataFrame(emis, SB.K)
        } 
        
        else if (class(emis) == "character") {
          switch (tools::file_ext(input,
            "csv" = private$readfromcsv(input, ...),
            "xlsx" = private$readFromExcel(input, ...)
          )) 
        }

        # else if ("unitFactor" %in% names(MoreParams)){
        #   private$UnitFactor <- MoreParams[["unitFactor"]]
        # }
      },
      
      emissions = function(scenario = NULL){
        if (is.NULL(private$emissions)) return (NULL)
        if (is.null(scenario)) 
          return (names(private$emissions))[!names(private$emissions) %in% c(
                                           "Scenario","VarName","Unit")]
      },
      
      CleanEmissions = function(value) { 
        if (missing(value)) {
          private$EmissionSource
        } else {
          
        }
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
      
      setEmissionFunction = function(app_input, kaas){
        browser()
        states <- colnames(kaas)
        
        if ("list" %in% class(app_input)) {
          # Check if the list was provided in the correct format
          if(!all(as.character(names(app_input)) %in% as.character(states))){
            stop("Abbreviations are incorrect")
          }
          for(i in app_input){
            if(!is.function(i)){
              stop("Not all elements in the list are functions")
            }
          }
          Emis <- app_input
          private$EmissionSource <- Emis
        } else if ("data.frame" %in% class(app_input)) {
          if(!(all(c("Abbr","Emis") %in% names(app_input)))){
            stop("Expected 'Abbr', 'Emis' and 'Timed' columns in dataframe")
          }
          if(!all(as.character(app_input$Abbr) %in% as.character(states))){
            stop("Abbreviations are incorrect")
          }
          Emis <- 
            app_input |> 
            group_by(Abbr) |> 
            summarise(n=n(),
                      EmisFun = list(
                        approxfun(
                          data.frame(Timed = c(0,Timed), 
                                     Emis=c(0,Emis)),
                          rule = 2)
                      )
            )
          
          funlist <- Emis$EmisFun
          names(funlist) <- Emis$Abbr
          
          private$EmissionSource <- funlist
          
        } else {
          stop("Expected a list of functions or dataframe with column 'Timed'")
        }
      },
      
      setEmissionDataFrame = function(emis_df, kaas) {
        #browser()
        if ("data.frame" %in% class(emis_df) && all(c("Abbr","Emis") %in% names(emis_df))) {
          #we need states - via solver from the core
          states <- colnames(kaas)
          # there can be multiple times in optional Timed column. 
          #if so, emis_df becomes a list of vectors
          if ("Timed" %in% names(emis_df)) {
            vEmis <- emis_df |>
              rename(var = Abbr) |>
              rename(time = Timed) |>
              rename(value = Emis) |>
              mutate(method = "add")
            
            private$EmissionSource <- vEmis
          } else { 
            #steady state
            vEmis <- rep(0.0, length.out = length(states))
            names(vEmis) <- states
            vEmis[match(emissions$Abbr, states)] <- emissions$Emis
            private$EmissionSource <- vEmis
            names(private$EmissionSource) <- states
            
          }
        }
      }
    )
  )