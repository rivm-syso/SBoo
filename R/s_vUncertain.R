#' @title vUncertain
#' @name vUncertain
#' @description sensitivity Mass 4 kaas,
#' #logged sd / intervals ??
#'  by changing kaas what is the change in Mass
#' @param ParentModule SBCore
#' @param vnamesDistSD  dataframe with columns vnames, distNames (see lhs package for possible distributions), secondPar
#' @param n samplesize 
#' @return States (i) (=mass)
#' @export
vUncertain = function(ParentModule, 
                      vnamesDistSD, n, tol=1e-30) { 
  
  TheCore <- ParentModule$myCore
  AllSBData <- TheCore$metaData()
  AllSBVars <- AllSBData$AttributeNames[!AllSBData$AttributeNames %in% names(TheCore$moduleList)]
  
  uniqvNames <- unique(vnamesDistSD$vnames)
  
  if (!all(uniqvNames %in% AllSBVars)) {
    stop(do.call(paste, c(list("Not all vnames found:"), 
                          as.list(uniqvNames[!uniqvNames %in% AllSBVars]))))
  }
  
  #fetch all variables, keep them as base
  baseVars <- lapply(uniqvNames, TheCore$fetchData)
  names(baseVars) <- uniqvNames
  
  # other preps; basic solution:
  SB.K = ParentModule$SB.k
  if (is.null(ParentModule$emissions))
    stop("No emissions in vUncertain")

  #emissions can be a list data.frame per State/t or just a vector
  if ("list" %in% class(ParentModule$emissions)) {
    firstEmis <- sapply(ParentModule$emissions, function(x) x$emis[1])
    basicSolv <- solve(SB.K, -firstEmis, tol = tol)
  } else {
    basicSolv <- solve(SB.K, -ParentModule$emissions, tol = tol)
  }

  unif01LHS <- lhs::optimumLHS(n=n, k=nrow(vnamesDistSD), maxSweeps=2, eps=.1, verbose=FALSE)
  #prep to save for analyses
  vnamesDistSD <- cbind(vnamesDistSD, t(unif01LHS))
  
  resultsAsList <- list()
  
  for (i in 1:n) {#loop hypercube
    
    #empty dataframe per uniq variable in vnamesDistSD
    Updated <- lapply(baseVars, function(x) {
      if ("data.frame" %in% class(x)) {
        x[F,]
      } else {
        unlist(unname(x))
      }
    }) 
    names(Updated) <- names(baseVars)

    for (vari in 1:nrow(vnamesDistSD)){ #loop vnamesDistSD rows 

      vname <- vnamesDistSD$vnames[vari]
      TakeDefault <- !("mean" %in% names(vnamesDistSD) && !is.na(vnamesDistSD$mean[vari]))
      #transform unif to scaling factor 
      scalingF <- switch (vnamesDistSD$distNames[vari],
                          "normal" = qnorm(p = unif01LHS[i, vari], mean = 1, sd = vnamesDistSD$secondPar[vari]),
                          "uniform" =  1 + vnamesDistSD$secondPar[vari] * (unif01LHS[i, vari] - 0.5)
      )
      
      if ("data.frame" %in% class(baseVars[[vname]])) {
        #apply to new row of data to Mutate
        Dees <- The3D[The3D %in% names(vnamesDistSD)]
        NeedDees <- Dees[!is.na(vnamesDistSD[vari,Dees])]
        if (TakeDefault) {
          newRowi <- which(sapply(NeedDees, function(aD){
            baseVars[[vname]][,aD] == vnamesDistSD[vari, aD]}))
          newRow <- baseVars[[vname]][newRowi, c(NeedDees, vname)]
          newRow[1,vname] <- newRow[1,vname] * scalingF
        } else {  # take mean from input
          newRow <- vnamesDistSD[vari, NeedDees]
          newRow[1,vname] <- vnamesDistSD$mean * scalingF 
        }
        Updated[[vname]] <- rbind(Updated[[vname]], newRow)
        
      } else { #just a number
        if(TakeDefault) {

          Updated[[vname]] <- scalingF * Updated[[vname]] 
        } else {
          Updated[[vname]] <- scalingF * vnamesDistSD$mean
        }
        #names(Updated[[vname]]) <- vname
      }
    }
    TheCore$UpdateVars(Updated)

    #update core and solve
    TheCore$UpdateDirty(uniqvNames)
    ParentModule$PrepKaasM()

    tryCatch(
      tryResult <- solve(ParentModule$SB.k, -ParentModule$emissions, tol = tol),
      error = function(e) {
        #collecting Updated and write to file
        #put column with varName in column "varName" and bind rows extending columns to all unique column names
        for (varname in names(Updated)) {
          if ("data.frame" %in% class(Updated[[varname]])) {
            Updated[[varname]]$varName <- varname
            names(Updated[[varname]])[names(Updated[[varname]]) == varname] <- "Waarde"
          } else { #just a listed value
            Updated[[varname]] <- data.frame(
              varName = varname,
              Waarde = unlist(Updated[[varname]])
            )
          }
        }
        write.csv(as.data.frame(do.call(bind_rows, Updated)), file = "testScripts/solvError.csv")
        stop(paste("Error in Solving; ", e))
    })
    
    resultsAsList[[as.character(i)]] <- tryResult
  }
  
  ParentModule$vnamesDistSD <- vnamesDistSD

  #reset variables to the originals
  for (i in 1: length(baseVars)){
    TheCore$SetConst(baseVars[[vname]])
  }
  TheCore$UpdateDirty(vnamesDistSD$vnames)
  
  #Save base as last, this reset eqMass
  resultsAsList[["base"]] <- basicSolv
  
  do.call(rbind, resultsAsList)
  
}
