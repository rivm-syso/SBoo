#'@title Resuspension rate constant of substances in sediment
#'@name k_Resuspenion
#'@description
#'Calculation of the resuspension rate constant from sediment to water
#'For molecular substanes. More documentation to be added.
#'@param to.SettlingVelocity SettlingVelocity, see function
#'@param VertDistance mixed depth sediment compartment #[m]
#'@param RhoCP Mineral density of sediment and soil #[kg/m3]
#'@param FRACs Volume fraction solids in sediment #[-]
#'@param to.mConcSusp Concentration of suspended matter in water #[kg/m3]
#'@return k_Resuspension Resuspension flow from sediment #[s-1]
#'@export

k_Resuspension <- function (VertDistance, #SettlVelocitywater
                            DynViscWaterStandard,
                            to.rhoMatrix,
                            to.NETsedrate,
                            to.RadCP, to.RhoCP, from.RhoCP, FRACs, to.SUSP, SpeciesName, ScaleName, to.SubCompartName, from.SubCompartName, Test) {

  if (SpeciesName == "Molecular"){
    if (as.character(Test) == "TRUE"){
      SettlingVelocitySPM <- 2.5/(24*3600)
    } else {
      # ScaleName
      SettlingVelocitySPM <-    f_SetVelWater(radius = to.RadCP,
                                              rhoParticle = to.RhoCP, rhoWater = to.rhoMatrix, DynViscWaterStandard)
    }
  } else {
    # ScaleName
  SettlingVelocitySPM <-    f_SetVelWater(radius = to.RadCP,
                                        rhoParticle = to.RhoCP, rhoWater = to.rhoMatrix, DynViscWaterStandard)
  }
  
  #Gross sedimentation rate from water [m/s]
  GROSSEDrate <- SettlingVelocitySPM*to.SUSP/(FRACs*from.RhoCP)    #[m.s-1] possibly < NETsedrate
  
  #Resuspension flow from sediment [m/s]; can't be < 0
  RESUSflow <- max(0, GROSSEDrate - to.NETsedrate) # for particulates this NETsedrate is not optimal!
  
  #Resuspension k to water [s-1]
  return(RESUSflow / VertDistance)
}



