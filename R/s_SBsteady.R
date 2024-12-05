#' @title SBsteady
#' @name SBsteady
#' @description run system of 1rst order k(i,j) and emissions,
#'  till steady state
#' @param ParentModule passing the SolverModule for access to it's public functions
#' @param tmax 
#' @return see: runsteady()
#' export
SBsteady = function(ParentModule, tmax=1e10) {
  
  SB.K = ParentModule$SB.k
  vEmis = ParentModule$emissions()
  
  m0 <- rep(0,nrow(SB.K))
  
  solved = rootSolve::runsteady(
    y = m0,
    times = c(0,tmax),
    func = ParentModule$SimpleBoxODE,
    parms = list(K = SB.K, e = vEmis)
  )
  
  if(!attr(solved, "steady")) {
    warning("steady state not reached")
  } else {
    inYear <- round(attr(solved, "time") / (60 * 60 * 24 * 365), 1) #365.25
    cat(paste("steady state reached at t =", inYear, "yr\n"))
  }
  
  names(solved)[names(solved) == "y"] <- "eqMass"
  names(solved)[names(solved) == "signal"] <- "emis"
  return(solved)
  
}