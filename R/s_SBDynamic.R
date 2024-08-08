#' @title SBsolve
#' @name SBDynamic
#' @description solve system of States for 1rst order k(i,j) and emissions
#' @param ParentModule SBcore object; to access SimpleBoxODE, states etc. 
#' too much and to diverse to put them in parameters and complicate the execute method
#' @param tmax time [s] for the simulation period
#' @param nTIMES number of timesteps
#' @param EmisAsPulse T is a pulse-type emission, false is not
#' @return see: ode()
  SBDynamic <- function(tmax = 1e10, nTIMES = 100) {
    
    SB.K = ParentModule$SB.k
    SBNames = colnames(SB.K)
    SB.m0 <- rep(0, length(SBNames))
    SBtime <- seq(0,tmax,length.out = nTIMES)
    funlist <- ParentModule$funlist
    
    out <- deSolve::ode(
      y = as.numeric(SB.m0),
      times = SBtime ,
      func = ODEapprox,
      parms = list(K = SB.K, SBNames=SBNames, funlist=funlist),
      rtol = 1e-10, atol = 1e-2)
    
    
    #if(as.character(class(deS)[1])!="data.frame") return (list(errorstate="error", deS))
    colnames(out)[1:length(SBNames)+1] <- SBNames
    colnames(out)[grep("signal",colnames(out))] <- paste("emis",SBNames,sep = "2")
    as.data.frame(out)
    
  
  deS
}