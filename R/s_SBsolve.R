#' @title SBsolve
#' @name SBsolve
#' @description solve system of States for 1rst order k(i,j) and emissions
#' @param tmax 
#' @return see: ode()
SBsolve = function(ParentModule, tmax = 1e10, nTIMES = 100) {

  SB.K = ParentModule$SB.k
  vEmis = ParentModule$emissions
  
  t <- seq(0,tmax,length.out = nTIMES)
  m0 <- rep(0,nrow(SB.K))
  
  deS <- deSolve::ode(
    y = m0,
    times = t,
    func = ParentModule$SimpleBoxODE,
    parms = list(K = SB.K, e = vEmis)
  )
  #if(as.character(class(deS)[1])!="data.frame") return (list(errorstate="error", deS))

  deS[,-1]
}