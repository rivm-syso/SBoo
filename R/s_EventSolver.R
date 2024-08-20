#' @title EventSolver
#' @name Eventsolver
#' @description solve system of States for 1rst order k(i,j) and emissions
#' @param ParentModule SBcore object; to access SimpleBoxODE, states etc. 
#' too much and to diverse to put them in parameters and complicate the execute method
#' @param tmax time [s] for the simulation period
#' @param nTIMES number of timesteps
#' @param EmisAsPulse T is a pulse-type emission, false is not
#' @return see: ode()
EventSolver = function(ParentModule, tmax = 1e10, nTIMES = 100) {
  
  SB.K = ParentModule$SB.k
  SBNames = colnames(SB.K)
  SB.m0 <- setNames(rep(0, length(SBNames)), SBNames)
  SBtime <- seq(0,tmax,length.out = nTIMES)
  vEmis <- ParentModule$emissions
  
  out <- deSolve::ode(
    y = SB.m0,
    times = SBtime ,
    func = ParentModule$EventODE,
    parms = list(K = SB.K, SBNames=SBNames),
    events = list(data = vEmis),
    rtol = 1e-11, atol = 1e-3)
  
  colnames(out)[1:length(SBNames)+1] <- SBNames
  as.data.frame(out) 
  
  attr(out, "SolverType") <- "EventSolver"
  
  return(out)
  
}