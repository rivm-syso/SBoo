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
                          kdeg.soil, kdeg.sed, Q.10,
                          kdeg.water, FRinw, BACTtest,BACTcomp, 
                          Ksw, Biodeg, CorgStandard, RHOsolid,
                          Matrix, SpeciesName) {
  
  f_kdegcalc(Q.10, Ksw, Biodeg, C.OHrad.n, Ea.OHrad, k0.OHrad,
                          CorgStandard, RHOsolid,
                          Matrix) 
  
  if (SpeciesName %in% c("Molecular")) {
    
    switch(Matrix,
           "air" =   {
             if (is.na(kdeg.air)){
               kdeg.air <-   f_kdegcalc(Q.10, Ksw, Biodeg, C.OHrad.n, Ea.OHrad, k0.OHrad,
                                        CorgStandard, RHOsolid,
                                        Matrix)
             }
             
             Tempfactor <- exp((Ea.OHrad/constants::syms$r)*((Temp-T25)/T25^2))
             FRingas * kdeg.air * (C.OHrad / C.OHrad.n) * Tempfactor
           },
           "soil" = {
             if (is.na(kdeg.soil)){
               kdeg.soil <-   f_kdegcalc(Q.10, Ksw, Biodeg, C.OHrad.n, Ea.OHrad, k0.OHrad,
                                        CorgStandard, RHOsolid,
                                        Matrix)
             }
             Tempfactor.wsds(Q.10,Temp,T25)*kdeg.soil
           },
           "sediment" = {
             if (is.na(kdeg.sed)){
               kdeg.sed <-   f_kdegcalc(Q.10, Ksw, Biodeg, C.OHrad.n, Ea.OHrad, k0.OHrad,
                                         CorgStandard, RHOsolid,
                                         Matrix)
             }
             Tempfactor.wsds(Q.10,Temp,T25)*kdeg.sed
           },
           "water" = {
             if (is.na(kdeg.water)){
               kdeg.water <-   f_kdegcalc(Q.10, Ksw, Biodeg, C.OHrad.n, Ea.OHrad, k0.OHrad,
                                         CorgStandard, RHOsolid,
                                         Matrix)
             }
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
             
             kdeg.sed # not corrected for temperature or other aspects
           },
           "water" = {
             kdeg.water # not corrected for temperature or other aspects
           }
    )
  }
  
}
