#' @title k_Degradation
#' @name k_Degradation
#' @description calculate k for degradation; 
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
k_Degradation <- function(FRingas, KdegDorC, C.OHrad.n, C.OHrad, 
                         Tempfactor,
                          FRinw, BACTtest,BACTcomp,
                          Matrix, SpeciesName, SubCompartName, ScaleName) {
  # exclusions of process:
  if (((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & (SubCompartName == "freshwatersediment" | 
                                                            SubCompartName == "lakesediment" |
                                                            SubCompartName == "lake" |
                                                            SubCompartName == "river" |
                                                            SubCompartName == "agriculturalsoil"|
                                                            SubCompartName == "othersoil")) | 
      (ScaleName %in% c("Regional", "Continental")) & (SubCompartName == "deepocean" )) {
    return(NA)
  }
    
  if (SpeciesName %in% c("Molecular")) {
    
    switch(Matrix,
           "air" =   {
             FRingas * KdegDorC * (C.OHrad / C.OHrad.n) * Tempfactor
           },
           "soil" = {
             Tempfactor*KdegDorC
           },
           "sediment" = {
             Tempfactor*KdegDorC
           },
           "water" = {
             Tempfactor*KdegDorC*(BACTcomp/BACTtest)*FRinw
           }
           
    )
  } else { # Particulate

    switch(Matrix,
           "air" =   {
             
             Tempfactor*KdegDorC # not corrected for temperature or other aspects
           },
           "soil" = {
             Tempfactor*KdegDorC # not corrected for temperature or other aspects
           },
           "sediment" = {
             
             Tempfactor*KdegDorC # not corrected for temperature or other aspects
           },
           "water" = {
             Tempfactor*KdegDorC # not corrected for temperature or other aspects
           }
    )
  }
  
}
