#'@title k_Deposition
#'@name k_Deposition
#'@param FRingas
#'@param WINDspeed
#'@param MW
#'@param MTCair
#'@param MTCother
#'@param twet
#'@param tdry
#'@param COLLECTeff
#'@param AEROSOLdeprate
#'@param Kaw
#'@param FRorig Fraction original species in water or porewater
#'@param FRaerw
#'@param FRaers 
#'@param SpeciesName Name of the species the function is applied to (now only implemented for molecular, particulate drydep to be added) [NA]
#'@param AREAFRAC Fraction of area compartment relative to the whole surface area of this scale (systemarea) [-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance Mixing depth air compartment #[m]
#'@param OtherkAir correction term for processes other than deposition affecting the air compartment
#'@return Transfer rate constant gaseous species air to soil or water #[s-1]
#'@export
k_Deposition <- function(FRingas, 
                         WINDspeed,
                         VertDistance, 
                         twet ,
                         tdry ,
                         COLLECTeff, # Can in SB5 be updated to (function from scavening in cloudwater stuf of SB4N!
                         AEROSOLdeprate , # should be transferred to data (or a function related to P species)
                         Kacompw,
                         FRorig,
                         SpeciesName,
                         OtherkAir,
                         landFRAC,
                         Kaers,
                         Kaerw,
                         FRACaer){ 
  
  if (SpeciesName %in% c("Molecular")) {
    
    HEIGHT.a <- VertDistance
    
    # NOT_deposition <- ((GASABS.a.w*(AREAFRAC.w0R+AREAFRAC.w1R+AREAFRAC.w2R)+
    #                       GASABS.a.s*(AREAFRAC.s1R+AREAFRAC.s2R+AREAFRAC.s3R))/HEIGHT.aR + 
    #                      KDEG.aR +
    #                      k.aR.aC) # problematic correction for other removal processes from air affection actual deposition.
    FRinaerw = fFRinaerw(Kaers, Kaerw, FRACaer)
      FRinaers = fFRinaers(Kaers, Kaerw, FRACaer)
      DRYDEPaerosol <- AEROSOLdeprate*(FRinaerw+FRinaers)
    AerosolWashout <- FRaers*(tdry+twet)/twet*COLLECTeff*RAINrate
    GasWashout <- FRingas*(tdry+twet)/twet*RAINrate/(Kacompw*FRorig)
    
    kdry <- DRYDEPaerosol/HEIGHT.a  + OtherkAir
    
    kwet <- (AerosolWashout+GasWashout)/HEIGHT.a + OtherkAir
    
    
    MeanRemAir <- ((1/kdry)*tdry/(tdry+twet)+(1/kwet)*twet/(tdry+twet)-
                     ((1/kwet-1/kdry)^2/(tdry+twet))*
                     (1-EXP(-kdry*tdry))*(1-EXP(-kwet*twet))/(1-EXP(-kdry*tdry-kwet*twet)))^-1
    
    
    
    
    MeanDep <- (MeanRemAir - OtherkAir) * landFRAC
    
    return( MeanDep ) # the gasabs here is for the two compartmenbts for which this function is run (air to soil/water)
    
  } else { # 
    return(NA)
  }
}
