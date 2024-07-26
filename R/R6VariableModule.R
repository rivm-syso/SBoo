#' @title Variable Module 
#' @description Module dealing with computing and initializing the variables
#' @import R6
#' @export
VariableModule <-
  R6::R6Class(
    "VariableModule",
    inherit = CalcGraphModule,
    public = list(
      execute = function(debugAt = NULL) {
        ret <- private$Execute(debugAt)
        if (any(is.na(ret))){
          return(data.frame(NA))
        } else { #an | for both cases would cause an error
          if (nrow(ret) == 0) {
            return(data.frame(NA))
          } else {
            return(ret)
          }
        }
      },
      
      initialize = function(TheCore, exeFunction, AggrBy, AggrFun){
        super$initialize(TheCore, exeFunction)
        private$aggrBy <- AggrBy
        private$aggrFUN <- AggrFun
      }
    ),
    active = list(
      aggr.by =  function(value) {
        if (missing(value)) {
          private$aggrBy
        } else {
          stop("`$aggr.by` is set at initialize()", call. = FALSE)
        }
      },
      aggr.FUN = function(value) {
          if (missing(value)) {
            private$aggrFUN
          } else {
            stop("`$aggr.FUN` is set at initialize()", call. = FALSE)
          }
        }
    ),
    
    private = list(
      aggrBy = NA,
      aggrFUN = NA,
      Execute = function(debugAt = NULL){ #debug can be list of (3D)names with values
        Fpars <- formalArgs(self$exeFunction)
        Fpars <- Fpars[Fpars!="..."] # leave out the obligatory?x param ...

        # Variable defining functions can have all. as pre.. indicating the need for the whole table
        all = as.logical(startsWith(unlist(Fpars), "all."))
        Fpars[all] <- substring(Fpars[all],5)
        
        #fetch input tables
        AllInput <- lapply(Fpars, private$MyCore$fetchData)
        #the constants separate
        InputIsDataFrame <- sapply(AllInput, function(x){
          "data.frame" %in% class(x)
        })
        AllConstants <- AllInput[!InputIsDataFrame]
        names(AllConstants) <- Fpars[!InputIsDataFrame]
        
        #and the all. type separate
        if (any(all)) {
          all.type <- AllInput[all]
          names(all.type) <- paste("all", Fpars[all], sep = ".")
        } else {
          all.type <- list() #it should exist
        }
        
        #the remainder, i.e. "normal" inputs
        AllInput <- AllInput[InputIsDataFrame & !all]
        
        #TODO remove this wizard mode
        #to debug the assembly of data 
        if ("assembly" %in% names(debugAt)) {
          browser()
        }
        #make dimension-scaffold
          #which dimensions?
        dims <- unique( unlist( lapply(AllInput, function(x){
          The3D[The3D %in% names(x)]
        })))
        To.3D <- paste("To.", The3D, sep = "")
        To.dims <- unique( unlist( sapply(AllInput, function(x){
          To.3D[To.3D %in% names(x)]
        })))
        if (length(To.dims) > 0 ) {
          stop ("Not possible to use 6 dims for a variable; process/flux only") #scaffold would be 2 big
        }
        #stopifnot(length(dims) > 0)
        if (length(dims) == 0) {# any possible exception will occur; all constants and/or all.-data.frames
          if (!is.null(debugAt)){
            debugonce(self$exeFunction)
          }
          ParsList <- as.list(c(AllConstants, all.type))
          ret <- data.frame(
            OneOnly = do.call(self$exeFunction, ParsList)
          ) 
          names(ret) <- private$MyName
          return(ret)

        } else {#regular case with SB variables having dimensions

          # merge the data frames (each like the SBcore$fetchtable result) 
          scaffold <- data.frame(unique(private$MyCore$states$asDataFrame[,dims,drop = F]))
          
          #TODO clean expand2scaffold, 
          expand2scaffold <- T
          if(expand2scaffold) {
            mergeScaffold <- function(x,y) {merge(x,y,all.x = T)}
            CalcTable <- Reduce(mergeScaffold, c(list(scaffold),AllInput))
          } else {
            CalcTable <- Reduce(merge, c(list(scaffold),AllInput))
            CalcTable <- CalcTable[complete.cases(CalcTable),]
          }
          
          DimColumns <- names(CalcTable) %in% The3D
          DimTable <- CalcTable[,names(CalcTable)[DimColumns],drop = F]
          
          #names(DimTable) <- names(CalcTable)[DimColumns] #lost if 1 dimension only, pffff
          CalcTable <- CalcTable[,names(CalcTable)[!DimColumns],drop = F]
          
          #prep debugnames for use in loop
          namesdebugAt <- names(debugAt)[names(debugAt) != "assembly"]
          if (length(debugAt) > 0 && names(debugAt) == "assembly") debugAt <- NULL # follow the apply case
          
          #needs debugging? If not: apply should be faster
          if (is.null(debugAt) && 
                !(length(CalcTable) != 1 | length(all.type) > 0)){
            
            if (length(CalcTable) == 1) {
              NewData <- do.call(lapply, 
                                 args = c(list(X = CalcTable[,1],
                                             FUN = self$exeFunction),
                                             AllConstants, all.type)
                                 )
                                #SIMPLIFY = F) #, USE.NAMES = USE.NAMES
            } else {
              argsVec <- lapply(1:length(CalcTable), function(x) CalcTable[,x])
              names(argsVec) <- names(CalcTable)
              if (length(AllConstants) > 0) {
                  NewData <- do.call("mapply",
                                     args = c(FUN = self$exeFunction, 
                                              argsVec, 
                                              MoreArgs = list(AllConstants),
                                              list(SIMPLIFY = F)))
              } else {
                  NewData <- do.call("mapply",
                                     args = c(FUN = self$exeFunction, 
                                              argsVec, 
                                              list(SIMPLIFY = F)))
              }
            }
          } else { #Call function for each row
            NewData <- lapply(1:nrow(CalcTable), function(i) {
              ParsList <- as.list(c(CalcTable[i,,drop = F], 
                                  AllConstants, 
                                  all.type))
            
              if (is.null(debugAt)) {
                ToDebug <- F
              } else {
                if(!all(namesdebugAt %in% names(ParsList))){
                  warning(paste("!all(namesdebugAt %in% names(ParsList)) in", private$MyName))
                }
                ToDebug <- T
                if (length(namesdebugAt) > 0 ){
                  for (j in 1:length(debugAt)){
                    if (ParsList[[names(debugAt)[j]]] != debugAt[j]){
                      ToDebug <- F
                      break
                    }
                  }
                }
              }
              #test if names of debugAT are in the ParsList
              if (ToDebug) debugonce(self$exeFunction)
            
              do.call(self$exeFunction, ParsList)
            })
          }
          CatchEmpty <- sapply(NewData, length) == 0
          ResultIsNA <- sapply(NewData, function(x) {
            is.na(x) | is.null(x)
          }) 
          CatchEmpty[!CatchEmpty] <- unlist(ResultIsNA)
          
          #remove the rows with NA results
          DimTable <- DimTable[!CatchEmpty,,drop=F]
          DimTable[[self$myName]] <-  unlist(NewData[!CatchEmpty])
          
          if (!is.na(private$aggrBy)) {
            #remove the aggrBy Dims 
            KeepNames <- private$aggrBy
            KeepNames <- KeepNames[KeepNames %in% names(DimTable)]
            
            if (is.na(private$aggrFUN)) {
              FUN <- sum
            } else {
              FUN <- match.fun(private$aggrFUN)
            }
            DimTable <- aggregate(DimTable[,self$myName], by = DimTable[,unlist(KeepNames),drop = F],
                                  FUN = FUN)
            #put correct names
            names(DimTable)[names(DimTable) == "x"] <- private$MyName
          }  
          
        }
        if (length(DimTable) == 0)
          browser()
        return(DimTable)

      }
    )
  )