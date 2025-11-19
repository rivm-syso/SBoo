#' @title Degradation rate constant measured or calculated
#' @name v_KdegDorC
#' @description calculate k for degradation for particulates and molecules.
#' if no specific value for kdeg is available it is scaled into categories 
#' as defined in Technical Guidance Document on risk assessment in support of
#' Commission Directive 93/67/EEC (European Commission, 2003)
#' @param kdeg degradation rate (as input) [s-1]
#' @param C.OHrad OH radical concentration specific to compartment, based on Wania & Daly (2002) [mol m-3]
#' @param C.OHrad.n general OH radical concentration, based on Wania & Daly (2002) [mol m-3]
#' @param k0.OHrad frequency factor of the OH radical reaction [m3 s-1] 
#' @param Ea.OHrad activation energy OH radical reaction [J mol-1]
#' @param T25 298K [K]
#' @param Q.10 rate increase factor per 10C [-]
#' @param KswDorC calculated soil water partitioning coefficient  [-]
#' @param BioDeg biodegradability test result [-]
#' @param CorgStandard Standard mass fraction organic carbon in soil/sediment [-]
#' @param rhoMatrix density of the matrix
#' @param Matrix compartment type considered 
#' @param SpeciesName species name considered
#' @return Degradation rate constant for molecular species
#' @export
KdegDorC <- function(kdeg, C.OHrad.n, k0.OHrad, Ea.OHrad, T25, 
                     Q.10, KswDorC, Biodeg, CorgStandard, rhoMatrix,
                     Matrix, SpeciesName,  Shortest_side, Intermediate_side, Longest_side,
                     Kssdr, Shape, CorFacSSA, RadS, deg_x, deg_tau, deg_y, 
                     deg_theta, deg_z, deg_eta, UVintensity, MICROBconc, Degrading_enzyme) {
  
  #Set default shortest, intermediate and longest side, in case it is not defined
  if ( is.na(Shortest_side) || is.null(Shortest_side) ) {
    Shortest_side <- RadS * 2
  }
  
  if ( is.na(Intermediate_side) || is.null(Intermediate_side) ) {
    Intermediate_side <- RadS * 2
  }
  
  if ( is.na(Longest_side) || is.null(Longest_side) ) {
    Longest_side <- RadS * 2
  }
  
  if (SpeciesName %in% c("Molecular")) {
    
    switch(Matrix,
           "sediment" = { 
             if(is.na(kdeg)){
               switch (Biodeg,
                       "r" = { a = # following table 5 in report 2015-0161
                         ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,30,
                                ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,300,
                                       ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,3000,
                                              ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,30000,NA))))}, # ready-biodegradable
                       "r-" = { a = 
                         ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,90,
                                ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,900,
                                       ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,9000,
                                              ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,90000,NA))))}, # ready-biodegradable (r-) substances failing the ten-day window
                       "i" = { a = 
                         ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,300,
                                ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,3000,
                                       ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,30000,
                                              ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,300000,NA))))}, # inherently biodegradable
                       "p" = { a = 
                         ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,300,
                                ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,3000,
                                       ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,30000,
                                              ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,300000,NA)))) }) # persistent
               return(0.1*Q.10^(13/10)*log(2)/a/(3600*24)) } else return(kdeg)
           },
           "soil" = {
             if(is.na(kdeg)){
               switch (Biodeg,
                       "r" = { a = # following table 5 in report 2015-0161
                         ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,30,
                                ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,300,
                                       ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,3000,
                                              ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,30000,NA))))}, # ready-biodegradable
                       "r-" = { a = 
                         ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,90,
                                ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,900,
                                       ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,9000,
                                              ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,90000,NA))))}, # ready-biodegradable (r-) substances failing the ten-day window
                       "i" = { a = 
                         ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,300,
                                ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,3000,
                                       ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,30000,
                                              ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,300000,NA))))}, # inherently biodegradable
                       "p" = { a = 
                         ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,300,
                                ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,3000,
                                       ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,30000,
                                              ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,300000,NA)))) }) # persistent
               return(Q.10^(13/10)*log(2)/a/(3600*24)) } else return(kdeg)
           },
           "water" = {
             if(is.na(kdeg)){
               switch (Biodeg,
                       "r" = { Q.10^(13/10)*log(2)/15/(3600*24) }, # ready-biodegradable
                       "r-" = { Q.10^(13/10)*log(2)/50/(3600*24) }, # ready-biodegradable (r-) substances failing the ten-day window
                       "i" = { Q.10^(13/10)*log(2)/150/(3600*24) }, # inherently biodegradable
                       "p" = { 1e-20 }) } else return(kdeg) # persistent
           },
           "air" = {
             if(is.na(kdeg)){
               return(C.OHrad.n * k0.OHrad * exp(-Ea.OHrad/(constants::syms$r*T25)))
             } else return(kdeg)
           })
  } else { 
    
    #Determine the shape factor, based on the approach of Maga et al. (2022) for surface degradation rate
    if (!is.na(Shape) && (Shape == "Cylindric - circular" | Shape == "Fiber")) {
      ShapeFac <- 3
      SAV <- 4/(Shortest_side*100)+2/(Longest_side*100)#Surface area to volume ratio - in cm-1
    } else if (!is.na(Shape) && Shape == "Film") {
      ShapeFac <- 2
      SAV <- 2/(Shortest_side*100)+2/(Intermediate_side*100)+2/(Longest_side*100)
    } else {
      ShapeFac <- 4 #If the shape is sphere, not given or irregular, use the Kssdr calculation for a sphere
      SAV <- 6/(Shortest_side*100)
    }
    
    #Calculate kdeg (s-1) using either of the following 3 approaches
    # Priority 1: empirical degradation parameters available (plasticFADE model)
    if (!any(is.na(c(deg_x, deg_tau, deg_y, deg_theta, deg_z, deg_eta,SAV)))) {
      kdeg <- deg_x * SAV^deg_tau * (deg_y * UVintensity^deg_theta + deg_z * MICROBconc^deg_eta) /(24*60*60) #plasticFADE results in a rate in d-1 to be converted in s-1
    # Priority 2: Kssdr available (Maga et al., 2022)
    } else if (!is.na(Kssdr)) {
      kdeg <- ShapeFac * Kssdr * CorFacSSA / (Shortest_side / 2)
      
    # Priority 3: default kdeg value
    } else {
      kdeg <- kdeg
    }
    
    #If polymer cannot be degraded in the absence of UV light (e.g. polyolefin) set kdeg to 0 in that compartment 
    if ((!is.na(Degrading_enzyme) && Degrading_enzyme == "FALSE") &&
        (!is.na(UVintensity) && UVintensity == 0)) {
      kdeg <- 0
    }
    
    switch(Matrix, #particulate
           "air" = kdeg,
           "soil" = kdeg,
           "sediment" = kdeg,
           "water" = kdeg,
           NA)
  }
  
}