#' @title k_Degradation
#' @name k_Degradation
#' @description calculate k for degradation in air; 
#' @param FRingas see
#' @param kdeg.air  kdeg.water kdeg.water
#' @param C.OHrad.n ?
#' @param C.OHrad ?
#' @param k0.OHrad ?
#' @param Ea.OHrad ?
#' @param Temp [K]
#' @param T25 = 298 K
#' @return Degradation rate constant for molecular species
#' @export
k_Degradation <- function(FRingas, kdeg.air, C.OHrad.n, C.OHrad, 
                          k0.OHrad, Ea.OHrad, Temp, T25, 
                          kdeg.soil, kdeg.sediment, Q.10,
                          kdeg.water, FRinw, BACTtest,BACTcomp, 
                          Matrix, SpeciesName) {
  
  if (SpeciesName %in% c("Molecular")) {
    
    switch(Matrix,
           "air" =   {
             Tempfactor <- exp((Ea.OHrad/constants::syms$r)*((Temp-T25)/T25^2))
             FRingas * kdeg.air * (C.OHrad / C.OHrad.n) * Tempfactor
           },
           "soil" = {
             Tempfactor.wsds(Q.10,Temp,T25)*kdeg.soil
           },
           "sediment" = {
             
             Tempfactor.wsds(Q.10,Temp,T25)*kdeg.sediment
           },
           "water" = {
             kdeg.water*Tempfactor.wsds(Q.10,Temp,T25)*(BACTcomp/BACTtest)*FRinw
           }
           
    )
  } else { # Particulate

    switch(Matrix,
           "air" =   {
             
             kdeg.air # not corrected for temperature or other aspects
           },
           "soil" = {
             kdeg.soil # not corrected for temperature or other aspects
           },
           "sediment" = {
             
             kdeg.sediment # not corrected for temperature or other aspects
           },
           "water" = {
             kdeg.water # not corrected for temperature or other aspects
           }
    )
  }
  
}
