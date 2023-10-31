#' @title k_Sedimentation
#' @name k_Sedimentation
#' @description Calculate the rate constant for sedimentation [s-1]
#' @param NaturalRho Density of natural particle
#' @param rhoMatrix Density of matrix in which particle is present
#' @param DynViscWaterStandard Dynamic viscosity of the fluid matrix
#' @param NaturalRad Radius of the natural particle
#' @param SettlingVelocity
#' @param CompartName function is defined for Water and Air; slightly different function
#' @param VertDistance compartment (mixed/mixk_Sedimentationing) DEPTH [m]
#' @return k_Sedimentation, the rate constant for sedimentation as first order process
#' @export
k_Sedimentation <- function(SettVellNat, FRinw, SettlingVelocity, DynViscWaterStandard,
                            SpeciesName, VertDistance, from.RhoCP, from.RadCP,
                            SubCompartName, ScaleName){
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & SubCompartName == "sea") {
    return(NA)
  } 
    # if (Matrix != "water") return (NA)
  
  if (SpeciesName == "Molecular") {
    SettVellNat <- f_SettlingVelocity(radius = from.RadCP,
                                      rhoParticle = from.RhoCP, rhoWater = 998, DynViscWaterStandard) 
    return(SettVellNat(1-FRinw) / VertDistance)
  } #else (return ends the function...)
  
  # SetVel <- f_SettlingVelocity (rad_species, rho_species, rhoMatrix, DynVisc,  
  #                               Matrix)
    
  SettlingVelocity/VertDistance
}