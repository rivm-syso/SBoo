#'@title k_Deposition
#'@name k_Deposition
#'@param FRingas
#'@param WINDspeed
#'@param MW
#'@param MTCair
#'@param MTCother
#'
#'@param MTC_2a
#'@param MTC_2w
#'@param MTC_2s

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
k_Deposition <- function(#MTC_2a,
                         #MTC_2s,
                         #MTC_2W,
                         #Kaw,
                         # from.ScaleName, 
                         # from.SubCompartName, 
                         # to.ScaleName,
                         # to.SubCompartName,
                         FRingas, 
                         VertDistance,
                         RAINrate,
                         twet ,
                         tdry ,
                         COLLECTeff, # Can in SB5 be updated to (function from scavening in cloudwater stuf of SB4N!
                         AEROSOLdeprate , # should be transferred to data (or a function related to P species)
                         Kacompw,
                         FRorig,
                         SpeciesName,
                         OtherkAir,
                         to.Area,
                         from.Area,
                         Kaers,
                         Kaerw,
                         FRinaerw,
                         FRinaers){ 
  
  if (SpeciesName %in% c("Molecular")) {
    
    DRYDEPaerosol <- AEROSOLdeprate*(FRinaerw+FRinaers)
    AerosolWashout <- FRinaers*(tdry+twet)/twet*COLLECTeff*RAINrate
    GasWashout <- FRingas*(tdry+twet)/twet*RAINrate/(Kacompw*FRorig) # (Kacompw * GRorig) == Aerosol collection efficiency
    
    kdry <- DRYDEPaerosol/VertDistance  + OtherkAir
    
    kwet <- (AerosolWashout+GasWashout)/VertDistance + OtherkAir
    
    
    MeanRemAir <- ((1/kdry)*tdry/(tdry+twet)+(1/kwet)*twet/(tdry+twet)-
                     ((1/kwet-1/kdry)^2/(tdry+twet))*
                     (1-exp(-kdry*tdry))*(1-exp(-kwet*twet))/(1-exp(-kdry*tdry-kwet*twet)))^-1
    
    MeanDep <- (MeanRemAir - OtherkAir) * (to.Area/from.Area)
    
    # The three lines below are work in progress. When this works, the degradation rate in R should be the same as in excel. Right now, the calculations of MeanDep and Degradation are different in R than in Excel.

    #GasAbs <- FRingas * (MTCair * MTCother/(MTCair * (Kaw * FRorig) + MTCother))
    
    #MeanDep <- (MeanRemAir - OtherkAir) 
    
    #Degradation <- (MeanDep + (GasAbs / VertDistance)) * (to.Area/from.Area)
    
    #return( Degradation )
    
    return( MeanDep ) # the gasabs here is for the two compartments for which this function is run (air to soil/water)
    
  } else { # 
    return(NA)
  }
}
