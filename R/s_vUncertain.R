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
vUncertain = function(ParentModule, vnamesDistSD, n, 
                      TargetScale = "Regional", TargetSubCompart = "river",
                      Targetspecies = "Molecular", tol=1e-30) { 
  
  TheCore <- ParentModule$myCore
  allves <- names(TheCore$moduleList)
  
  if (!all(vnamesDistSD$vnames %in% allves)) {
    stop(do.call(paste, c(list("Not all vnames found:"), 
                          as.list(vnamesDistSD$vnames[!vnamesDistSD$vnames %in% allves]))))
  }
  
  #fetch all variables, keep them as base
  baseVars <- lapply(vnamesDistSD$vnames, TheCore$fetchData)
  names(baseVars) <- vnamesDistSD$vnames
  
  # other preps; basic solution:
  SB.K = ParentModule$SB.k
  if (is.null(ParentModule$emissions))
    stop("No emissions in vUncertain")
  basicSolv <- solve(SB.K, -ParentModule$emissions, tol = tol)

  unif01LHS <- lhs::optimumLHS(n=n, k=nrow(vnamesDistSD), maxSweeps=2, eps=.1, verbose=FALSE)
  #prep to save for analyses
  vnamesDistSD <- cbind(vnamesDistSD, t(unif01LHS))
  
  aslist <- list()
  
  for (i in 1:n) {
    for (vari in 1:nrow(vnamesDistSD)){
      vname <- vnamesDistSD$vnames[vari]
      Updated <- baseVars[[vname]]
      
      #transform unif to scaling factor 
      scalingF <- switch (vnamesDistSD$distNames[vari],
        "normal" = qnorm(p = unif01LHS[i, vari], mean = 1, sd = vnamesDistSD$secondPar[vari]),
        "uniform" =  1 + vnamesDistSD$secondPar[vari] * (unif01LHS[i, vari] - 0.5)
      )
      vnamesDistSD[vnamesDistSD$vnames == vname, as.character(i)] <- scalingF
      Updated[,vname] <- scalingF * Updated[,vname]
      
      asParam <- list(Updated)
      names(asParam) <- vname
      do.call(TheCore$SetConst, asParam)
      
    }
    #update core and solve
    TheCore$UpdateDirty(vnamesDistSD$vnames)
    ParentModule$PrepKaasM()

    aslist[[as.character(i)]] <- solve(ParentModule$SB.k, -ParentModule$emissions, tol = tol)
  }
  
  ParentModule$vnamesDistSD <- vnamesDistSD

  #reset variables to the originals
  for (i in 1: length(baseVars)){
    TheCore$SetConst(baseVars[[vname]])
  }
  TheCore$UpdateDirty(vnamesDistSD$vnames)
  
  #Save base as last, this reset eqMass
  aslist[["base"]] <- basicSolv
  
  do.call(rbind, aslist)
  
}
