#' @title states
#' @description object in both SBcore as in SolverModule (slimmed down to existing kaas). 
#' @import R6
#' @export
SBstates <- R6::R6Class("SBstates",
  public = list(
    #' @description init, based on 
    #' @param StatesAsDataFrame dataframe of states with columns The3D Abbr and counter i
    initialize = function(ScaleSheet, SubCompartSheet, SpeciesSheet) {
      
      if (!all("data.frame" %in% sapply(list(ScaleSheet, SubCompartSheet, SpeciesSheet), class))) {
        stop("all inputs should have class data.frame")
      }
      self$buildStates(ScaleSheet, SubCompartSheet, SpeciesSheet)
    },
    
    #' @description Build the States table from the three dimension sheets
    #' @param ScaleSheet data.frame with at least columns: Scale, AbbrS, Default
    #' @param SubCompartSheet data.frame with at least columns: SubCompart, AbbrC, Default
    #' @param SpeciesSheet data.frame with at least columns: Species, AbbrP, Default
    #' @return Invisibly, the constructed States data.frame; also stores it in private$AsDataFrame
    buildStates = function(ScaleSheet, SubCompartSheet, SpeciesSheet) {
      # basic input checks
      required_scale_cols     <- c("Scale", "AbbrS", "Default")
      required_subcomp_cols   <- c("SubCompart", "AbbrC", "Default")
      required_species_cols   <- c("Species", "AbbrP", "Default")
      
      if (!all(required_scale_cols %in% names(ScaleSheet))) {
        stop("ScaleSheet must contain columns: ",
             paste(required_scale_cols, collapse = ", "))
      }
      if (!all(required_subcomp_cols %in% names(SubCompartSheet))) {
        stop("SubCompartSheet must contain columns: ",
             paste(required_subcomp_cols, collapse = ", "))
      }
      if (!all(required_species_cols %in% names(SpeciesSheet))) {
        stop("SpeciesSheet must contain columns: ",
             paste(required_species_cols, collapse = ", "))
      }
      
      ScaleSheet$Default       <- to_logical(ScaleSheet$Default)
      SubCompartSheet$Default  <- to_logical(SubCompartSheet$Default)
      SpeciesSheet$Default     <- to_logical(SpeciesSheet$Default)
      
      # select default entries per dimension
      default_scales      <- ScaleSheet$Scale[ScaleSheet$Default]
      default_subcomparts <- SubCompartSheet$SubCompart[SubCompartSheet$Default]
      default_species     <- SpeciesSheet$Species[SpeciesSheet$Default]
      
      if (length(default_scales) == 0L ||
          length(default_subcomparts) == 0L ||
          length(default_species) == 0L) {
        stop("No defaults found in one or more sheets (Scale/SubCompart/Species).")
      }
      
      # create all permutations (Cartesian product)
      States <- expand.grid(
        Scale      = default_scales,
        SubCompart = default_subcomparts,
        Species    = default_species,
        stringsAsFactors = FALSE
      )
      
      # attach abbreviations
      # match each dimension to its abbreviation
      abbrC <- SubCompartSheet$AbbrC[
        match(States$SubCompart, SubCompartSheet$SubCompart)
      ]
      abbrS <- ScaleSheet$AbbrS[
        match(States$Scale, ScaleSheet$Scale)
      ]
      abbrP <- SpeciesSheet$AbbrP[
        match(States$Species, SpeciesSheet$Species)
      ]
      
      if (any(is.na(abbrC))) {
        stop("Missing SubCompart abbreviations (AbbrC) for some states.")
      }
      if (any(is.na(abbrS))) {
        stop("Missing Scale abbreviations (AbbrS) for some states.")
      }
      if (any(is.na(abbrP))) {
        stop("Missing Species abbreviations (AbbrP) for some states.")
      }
      
      States$Abbr <- paste0(abbrC, abbrS, abbrP)
      
      # ensure column order: Abbr, Scale, SubCompart, Species
      States <- States[, c("Abbr", The3D), drop = FALSE]
      
      # store in private$AsDataFrame
      private$AsDataFrame <- States
      
      invisible(States)
    },
    
    #' @description Filter the States table using a set of rule rows
    #' @param rules data.frame with columns: Scale, SubCompart, Species, Operator, Value, Keep
    #' @return Invisibly, the filtered States data.frame; also updates private$AsDataFrame
    filterStates = function(rules) {
      
      if (is.null(private$AsDataFrame)) {
        stop("No States table available; call buildStates() first.")
      }
      
      States <- private$AsDataFrame
      
      required_rule_cols <- c("Scale", "SubCompart", "Species", "Operator", "Keep")
      if (!all(required_rule_cols %in% names(rules))) {
        stop("rules must contain columns: ", paste(required_rule_cols, collapse = ", "))
      }
      
      # start with all rows allowed
      keep_mask <- rep(TRUE, nrow(States))
      
      # helper to parse Value for %in%
      parse_value <- function(val, op = "%in%") {
        if (op == "%in%") {
          # comma-separated list -> character vector
          strsplit(as.character(val), "\\s*,\\s*")[[1]]
        } else {
          as.character(val)
        }
      }
      
      # apply each rule in turn
      for (i in seq_len(nrow(rules))) {
        
        rule <- rules[i, , drop = FALSE]
        
        op   <- rule$Operator
        keep <- isTRUE(rule$Keep)
        
        if (!op %in% c("==", "!=", "%in%")) {
          stop("Unsupported Operator in rules: ", op)
        }
        
        # build match mask for this rule: combine all non-NA dimensions with AND
        rule_mask <- rep(TRUE, nrow(States))
        
        # Scale condition
        if (!is.na(rule$Scale) && rule$Scale != "") {
          if (op == "==") {
            rule_mask <- rule_mask & (States$Scale == rule$Scale)
          } else if (op == "!=") {
            rule_mask <- rule_mask & (States$Scale != rule$Scale)
          } else if (op == "%in%") {
            rule_mask <- rule_mask & (States$Scale %in% parse_value(rule$Scale))
          }
        }
        
        # SubCompart condition
        if (!is.na(rule$SubCompart) && rule$SubCompart != "") {
          if (op == "==") {
            rule_mask <- rule_mask & (States$SubCompart == rule$SubCompart)
          } else if (op == "!=") {
            rule_mask <- rule_mask & (States$SubCompart != rule$SubCompart)
          } else if (op == "%in%") {
            rule_mask <- rule_mask & (States$SubCompart %in% parse_value(rule$SubCompart))
          }
        }
        
        # Species condition
        if (!is.na(rule$Species) && rule$Species != "") {
          if (op == "==") {
            rule_mask <- rule_mask & (States$Species == rule$Species)
          } else if (op == "!=") {
            rule_mask <- rule_mask & (States$Species != rule$Species)
          } else if (op == "%in%") {
            rule_mask <- rule_mask & (States$Species %in% parse_value(rule$Species))
          }
        }
        
        # combine with global keep_mask
        if (keep) {
          # keep only rows that either already kept AND match this rule
          keep_mask <- keep_mask & rule_mask
        } else {
          # drop rows that match this rule
          keep_mask <- keep_mask & !rule_mask
        }
      }
      
      States_filtered <- States[keep_mask, , drop = FALSE]
      
      private$AsDataFrame <- States_filtered
      
      invisible(States_filtered)
    },
    
    #' @description Filter the States table using mode-specific rules
    #' @param mode scalar (character/numeric) indicating the mode to apply
    #' @param mode_rules data.frame with columns: mode, Scale, SubCompart, Species, Keep
    #' @return Invisibly, the filtered States data.frame; also updates private$AsDataFrame
    filterStatesByMode = function(mode, mode_rules) {
      
      if (is.null(private$AsDataFrame)) {
        stop("No States table available; call buildStates() first.")
      }
      
      if (missing(mode)) {
        stop("filterStatesByMode() requires a 'mode' argument")
      }
      
      States <- private$AsDataFrame
      
      required_rule_cols <- c("mode", "Scale", "SubCompart", "Species", "Keep")
      if (!all(required_rule_cols %in% names(mode_rules))) {
        stop("mode_rules must contain columns: ",
             paste(required_rule_cols, collapse = ", "))
      }
      
      # subset rules to the requested mode
      rules <- mode_rules[mode_rules$mode == mode, , drop = FALSE]
      if (nrow(rules) == 0L) {
        # nothing to do; return current States unchanged
        return(invisible(States))
      }
      
      # start with all rows allowed
      keep_mask <- rep(TRUE, nrow(States))
      
      # apply each rule in turn
      for (i in seq_len(nrow(rules))) {
        
        rule <- rules[i, , drop = FALSE]
        keep <- isTRUE(rule$Keep)
        
        # build match mask for this rule: combine all non-NA/non-empty dimensions with AND
        rule_mask <- rep(TRUE, nrow(States))
        
        # Scale condition
        if (!is.na(rule$Scale) && rule$Scale != "") {
          rule_mask <- rule_mask & (States$Scale == rule$Scale)
        }
        
        # SubCompart condition
        if (!is.na(rule$SubCompart) && rule$SubCompart != "") {
          rule_mask <- rule_mask & (States$SubCompart == rule$SubCompart)
        }
        
        # Species condition
        if (!is.na(rule$Species) && rule$Species != "") {
          rule_mask <- rule_mask & (States$Species == rule$Species)
        }
        
        # combine with global keep_mask
        if (keep) {
          # keep only rows that already kept AND match this rule
          keep_mask <- keep_mask & rule_mask
        } else {
          # drop rows that match this rule
          keep_mask <- keep_mask & !rule_mask
        }
      }
      
      States_filtered <- States[keep_mask, , drop = FALSE]
      
      private$AsDataFrame <- States_filtered
      
      invisible(States_filtered)
    },
    
    #' @description map vector abbr for scale species subcompart to index in State
    #' @param abbr Oldschool abbreviation (character vector)
    #' @return integer vector of indices in private$AsDataFrame$Abbr (NA if not found)
    findState = function(abbr) {
      
      if (missing(abbr)) {
        stop("findState() requires an 'abbr' argument")
      }
      if (!is.character(abbr)) {
        stop("findState() expects 'abbr' to be a character vector")
      }
      
      # position of Scale (pattern must be one of R, C, A, M, T)
      m <- regexpr("[RCAMT]", abbr)
      
      # check for missing pattern
      if (any(m == -1)) {
        bad <- abbr[m == -1]
        warning(
          "Some abbreviations do not contain a scale code [RCAMT]: ",
          paste(bad, collapse = ", ")
        )
      }
      
      # helper to be safe with substr when m == -1
      safe_substr <- function(x, start, stop) {
        # if start or stop <= 0 or start > nchar(x), return ""
        if (length(x) == 0) return(character(0))
        n <- nchar(x)
        start[is.na(start)] <- 1
        stop[is.na(stop)] <- n
        out <- character(length(x))
        for (i in seq_along(x)) {
          if (start[i] < 1 || stop[i] < start[i] || start[i] > n[i]) {
            out[i] <- ""
          } else {
            out[i] <- substr(x[i], start[i], stop[i])
          }
        }
        out
      }
      
      # extract parts
      ScaleAbbr <- safe_substr(abbr, m, m)
      SUbCompartAbbr <- safe_substr(abbr, 1, m - 1)
      Spec <- safe_substr(abbr, m + 1, m + 1)
      
      # if all Spec equal "" it's data from the Molecular (pre Nano) version -> set to U
      empty_spec <- (Spec == "" | m == -1)
      if (all(empty_spec)) {
        Spec <- rep("U", length(Spec))
      } else {
        Spec[empty_spec] <- "U"  # for individual problematic entries
      }
      
      # G=Gas, D=Dissolved -> U=unbound
      Spec[Spec %in% c("D", "G")] <- "U"
      
      # harmonise subcompartment abbreviations
      SUbCompartAbbr[SUbCompartAbbr == "s"]  <- "s3"
      SUbCompartAbbr[SUbCompartAbbr == "sd"] <- "sd2"
      
      # build lookup strings
      lookup_abbr <- paste(SUbCompartAbbr, ScaleAbbr, Spec, sep = "")
      
      # match against the states table
      idx <- match(lookup_abbr, private$AsDataFrame$Abbr)
      
      # warn on any non-matched abbreviations
      if (any(is.na(idx))) {
        bad <- abbr[is.na(idx)]
        warning(
          "Some abbreviations could not be matched to states: ",
          paste(bad, collapse = ", ")
        )
      }
      
      idx
    },
    
    #' @description Change any of the 3D columns into ordered factors and sort a data.frame accordingly
    #' @param aDFwithD A data.frame to sort; must contain at least one of The3D or column 'Abbr'
    #' @return Sorted data.frame, with The3D columns converted to ordered factors using Dlevels
    sortFactors = function(aDFwithD) {
      
      if (!is.data.frame(aDFwithD)) {
        stop("sortFactors() expects a data.frame")
      }
      
      # detect which 3D columns are present
      Dcolumns <- The3D[The3D %in% names(aDFwithD)]
      
      # if no 3D columns, but Abbr is present, merge to get them
      if (length(Dcolumns) == 0 && "Abbr" %in% names(aDFwithD)) {
        aDFwithD <- merge(aDFwithD, private$AsDataFrame, by = "Abbr")
        Dcolumns <- The3D[The3D %in% names(aDFwithD)]
      }
      
      if (length(Dcolumns) == 0) {
        stop(
          "data.frame offered to sortFactors() must contain at least one of ",
          paste(The3D, collapse = ", "), " or column 'Abbr'"
        )
      }
      
      # replace the dimension columns with ordered factors using Dlevels
      for (theD in Dcolumns) { # The3D = c(Scale, SubCompart, Species)
        theName <- paste0(theD, "Name")
        levs <- self$Dlevels[[theD]]
        
        if (is.null(levs)) {
          stop("No levels found for dimension '", theD, "'. Is myCore set correctly?")
        }
        if (!all(c(theD, theName) %in% names(levs))) {
          stop("Levels table for '", theD,
               "' must contain columns '", theD, "' and '", theName, "'.")
        }
        
        aDFwithD[[theD]] <- factor(
          aDFwithD[[theD]],
          levels  = levs[[theD]],
          labels  = levs[[theName]],
          ordered = TRUE
        )
      }
      
      # sort by the 3D columns present (which are now ordered factors)
      ord <- do.call(order, unname(aDFwithD[Dcolumns]))
      aDFwithD[ord, , drop = FALSE]
    },
    
    #' @description find an element in any dimension, or the name of a dimension NOT FINISHED
    #' @param aDimension a dimension or an element therein
    #' @return list(the found dimension = the member) or list(the found dimension = dimension). 
    #' The found dimension is 1 of the 3dim
    findDim = function (aDimension) {
      browser()
      InDims <- sapply(The3D, function(aDim) {
        match(aDimension, private$AsDataFrame[,aDim])
      })
      if (sum(InDims>0)){
        
      }
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
    },
    myCore = function(value) {
      if (missing(value)) {
        return(private$MyCore)
      } else {
        if ("SBcore" %in% class(value)) {
          private$MyCore = value
        } else {
          stop("Core object expected, but not provided")
        }
      }
    },
    Dlevels = function(value) {
      if (missing(value)) {
        if (is.null(private$dlevels)) {
          #make/save levels to factors (sorting, in ggplot)
          ScaleOrder <- private$MyCore$fetchData("ScaleSheet")
          SubCompartOrder <- private$MyCore$fetchData("SubCompartSheet")
          SpeciesOrder <- private$MyCore$fetchData("SpeciesSheet")
          private$dlevels <- list(
            Scale = ScaleOrder[order(ScaleOrder$ScaleOrder), c("Scale", "ScaleName")],
            SubCompart = SubCompartOrder[order(SubCompartOrder$SubCompartOrder), c("SubCompart", "SubCompartName")],
            Species = SpeciesOrder[order(SpeciesOrder$SpeciesOrder), c("Species", "SpeciesName")]
          )
        }
        return(private$dlevels)
      } else {
        warning("Dlevels cannot be set; late initialized")
      }
    }
    
  ),
  private = list(
    AsDataFrame = NULL,
    MyCore = NULL,
    dlevels = NULL
  )
)