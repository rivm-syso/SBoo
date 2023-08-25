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
  vEmis = ParentModule$emissions
  
  m0 <- rep(0,nrow(SB.K))
  
  rootSolve::runsteady(
    y = m0,
    times = c(0,tmax),
    func = ParentModule$SimpleBoxODE,
    parms = list(K = SB.K, e = vEmis)
  )
}