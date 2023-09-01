#' @title states
#' @description object in both SBcore as in SolverModule (slimmed down to existing kaas). 
#' @import R6
#' @export
SBstates <- R6::R6Class("SBstates",
  public = list(
    #' @description init, based on 
    #' @param StatesAsDataFrame dataframe of states with columns The3D Abbr and counter i
    initialize = function(StatesAsDataFrame) {
      if (!all(The3D %in% names(StatesAsDataFrame))) {
        stop(paste("states should contain all of", The3D))
      }
      if (!all(c("Abbr") %in% names(StatesAsDataFrame))) {
        stop("states should contain Abbr")
      }
      if (!"i" %in% names(StatesAsDataFrame)) StatesAsDataFrame$i <- 1:nrow(StatesAsDataFrame)
      private$AsDataFrame <- StatesAsDataFrame[,c("i", "Abbr", The3D)]
    },
    
    #' @description return match (vector of i) in this states. i is the very original order of states
    #' @param vi with i in order of vi 
    matchi = function(vi) {
      if (is.null(private$AsDataFrame )) stop("states not properly initialized")
      match(vi, private$AsDataFrame$i)
    },
    
    #' @description map vector abbr for scale species subcompart to index in State
    #' @param abbr Oldschool abbreviation
    #' @return vector of indices
    findState = function (abbr) {
      #split at the centre == Scale
      m<-regexpr("[RCAMT]",abbr)
      ScaleAbbr <- substr(abbr,m,m)
      SUbCompartAbbr <- substr(abbr,1,m-1)
      #identical functioning; thus rename: s = soil -> s1 naturalsoil, sd = oceansediment -> sd2
      #lengthabbC <- sapply(CompartAbbr, nchar) #needed to detect missing Species
      SUbCompartAbbr[SUbCompartAbbr=="s"] <- "s3"
      SUbCompartAbbr[SUbCompartAbbr=="sd"] <- "sd2"
      Spec <- substr(abbr,m+1,m+1)
      #if all Spec equal "" it's data from the Molecular (pre Nano) version
      if (all(Spec == "")){ #put all to Unbound
        Spec <- rep("U", length(Spec))
      }
      #G=Gas, D=Dissolved -> U=unbound,
      Spec[Spec %in% c("D","G")] <- "U"
      PasteAndMatch <- function (abbrnr) {
        match(paste(SUbCompartAbbr[abbrnr],
                    ScaleAbbr[abbrnr],
                    Spec[abbrnr],sep=""), private$AsDataFrame$Abbr)
      }
      sapply(1:length(abbr),FUN = PasteAndMatch)
    }
    
  ),
  active = list(
    #' @field asDataFrame convienent returning states in a data.frame
    asDataFrame = function(value) {
      if (missing(value)) {
        private$AsDataFrame
      } else {
        stop("`$states` are set by new()", call. = FALSE)
      }
    },
    #' @field nStates just 4 convenience
    nStates = function(value) {
      if (missing(value)) {
        if (is.null(private$AsDataFrame)) return(0)
          nrow(private$AsDataFrame)
      } else {
        stop("property nStates is read-only", call. = FALSE)
      }
    }
    
    
  ),
  private = list(
    AsDataFrame = NULL
  )
)