#'@title Resuspension rate constant of substances in sediment
#'@name k_Resuspenion
#'@description
#'Calculation of the resuspension rate constant from sediment to water
#'For molecular substanes. More documentation to be added.
#'@param to.SettlingVelocity SettlingVelocity, see function
#'@param VertDistance mixed depth water sediment compartment #[m]
#'@param RHOsolid Mineral density of sediment and soil #[kg/m3]
#'@param FRACs Volume fraction solids in sediment #[-]
#'@param to.mConcSusp Concentration of suspended matter in water #[kg/m3]
#'@return k_Resuspension Resuspension flow from sediment #[s-1]
#'@export

k_Resuspension <- function (VertDistance, #SettlVelocitywater
                            rad_species, rho_species, to.rhoMatrix, DynViscWaterStandard,
                            to.Matrix,to.NETsedrate,
                            RHOsolid, FRACs, to.SUSP, SpeciesName) {
  
  SettlingVelocity2 <- switch (SpeciesName,
                              # "Molecular" = SettlVelocitywater, In the molecular version
                              "Molecular" =     f_SetVelWater(radius = from.RadCP,
                                                              rhoParticle = from.RhoCP, rhoWater = 998, DynViscWaterStandard) , # In the nano version
                              #else
                              SettlingVelocity
  )
  #Gross sedimentation rate from water [m/s]
  GROSSEDrate <- SettlingVelocity2*to.SUSP/(FRACs*RHOsolid)    #[m.s-1] possibly < NETsedrate
  
  #Resuspension flow from sediment [m/s]; can't be < 0
  RESUSflow <- max(0, GROSSEDrate - to.NETsedrate) # for particulates this NETsedrate is not optimal!
  
  #Resuspension k to water [s-1]
  return(RESUSflow / VertDistance)
}



