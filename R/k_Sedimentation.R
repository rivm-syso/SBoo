#' @title k_Sedimentation
#' @name k_Sedimentation
#' @description Calculate the rate constant for sedimentation [s-1]
#' @param FRinw Fraction chemical dissolved in water [-]
#' @param VertDistance Compartment depth/height [m]
#' @param DynViscWaterStandard Dynamic viscosity of the fluid matrix
#' @param RadCP Radius of the Coarse natural particle [m]
#' @param RhoCP Density of the Coarse natural particle [m]
#' @param SettlingVelocity Settling velocity of particulate species [m.s-1]
#' @param SubCompartName Name of relevant subcompartment for which k_Sedimentation is being calculated
#' @param ScaleName Name of relevant scale for which k_Sedimentation is being calculated
#' @param SpeciesName Name of relevant species (Molecular or particulate) for which k_Sedimentation is being calculated
#' @return k_Sedimentation, the rate constant for sedimentation as first order process
#' @export
k_Sedimentation <- function(FRinw, SettlingVelocity, DynViscWaterStandard,
                            VertDistance, from.RhoCP, from.RadCP,
                            SpeciesName, SubCompartName, to.SubCompartName, ScaleName){
  # if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & SubCompartName == "sea") {
  #   return(NA)
  # }
  if ((ScaleName %in% c("Regional", "Continental")) & to.SubCompartName == "deepocean") {
    return(NA)
  }
    # if (Matrix != "water") return (NA)
  
  if (SpeciesName == "Molecular") {
    SetlingVelocityCP <- f_SetVelWater(radius = from.RadCP,
                                      rhoParticle = from.RhoCP, rhoWater = 998, DynViscWaterStandard) 
    return(SetlingVelocityCP*(1-FRinw) / VertDistance)
  } else SettlingVelocity/VertDistance
}