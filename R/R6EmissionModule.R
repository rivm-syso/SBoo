#' @title Emission module
#' @description Function to get emissions from excel/csv files
#' @import R6
#' @export
EmissionModule <-
  R6::R6Class(
    "EmissionModule",
    public = list(
      initialize = function(TheCore, ...) {      
        private$MyCore <- TheCore
        private$MoreParams <- list(...)
        public$execute()
      },
      
      execute = function() {
        if (length(private$MoreParams) == 1 && class(private$MoreParams[[1]]) == "character") {
          #if it's a csv or xlsx file, we try reading it
          fn <- private$MoreParams[[1]]
          if (file.exists(fn)) {
            ext <- private$getFileNameExtension(fn)
            if (ext == "csv") {
              df <- read.csv(fn)
            } else if (ext == "xlsx") {
              df <- private$readFromExcel(fn)
            } else stop(paste("Cannot read", fn))
            
          } else stop(paste("file not found", fn))
        }
        return(private$Emissions)
      },
      
      emissions = function(scenario = NULL){
        if (is.NULL(private$emissions)) return (NULL)
        if (is.null(scenario)) 
          return (names(private$emissions))[!names(private$emissions) %in% c(
                                           "Scenario","VarName","Unit")]
      }
      
    ),

    private = list(
      MyCore = NULL,
      MoreParams = NULL,
      Emissions = NULL,
      BaseEmissions = NULL, # for statistical assessment
      getFileNameExtension = function (fn) { #grabfromtheweb thx pisca46
        # remove a path
        splitted    <- strsplit(x=fn, split='/')[[1]]   
        # or use .Platform$file.sep in stead of '/'
        fn          <- splitted [length(splitted)]
        ext         <- ''
        splitted    <- strsplit(x=fn, split='\\.')[[1]]
        l           <-length (splitted)
        if (l > 1 && sum(splitted[1:(l-1)] != ''))  ext <-splitted [l] 
        # the extention must be the suffix of a non-empty name    
        ext
      },
      readFromExcel = function(fn) {
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
      }
    )
  )