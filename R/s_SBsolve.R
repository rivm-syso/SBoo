#' @title SBsolve
#' @name SBsolve
#' @description solve system of States for 1rst order k(i,j) and emissions
#' @param ParentModule SBcore object; to access SimpleBoxODE, states etc. 
#' too much and to diverse to put them in parameters and complicate the execute method
#' @param tmax time [s] for the simulation period
#' @param nTIMES number of timesteps
#' @param EmisAsPulse T is a pulse-type emission, false is not
#' @return see: ode()
SBsolve = function(ParentModule, tmax = 1e10, nTIMES = 100, EmisAsPulse = F) {

  SB.K = ParentModule$SB.k
  
  SBtime <- seq(0,tmax,length.out = nTIMES)
  ParentModule$SBtime.tvars <- list(SBtime = SBtime)
  
  if (EmisAsPulse) {
    m0 = ParentModule$emissions
    vEmis =  rep(0,nrow(SB.K))
  } else { 
    #The default behaviour
    m0 <- rep(0,nrow(SB.K))
    vEmis= ParentModule$emissions
  }
  
  deS <- deSolve::ode(
    y = m0,
    times = SBtime,
    func = ParentModule$SimpleBoxODE,
    parms = list(K = SB.K, e = vEmis)
  )
  #if(as.character(class(deS)[1])!="data.frame") return (list(errorstate="error", deS))

  deS[,-1]
}