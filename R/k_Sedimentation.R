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
#' @param Test determines if SB4-Excel approach is taken or enhanced method from R version [boolean]
#' @return k_Sedimentation, the rate constant for sedimentation as first order process
#' @export
k_Sedimentation <- function(FRinw, SettlingVelocity, DynViscWaterStandard, rhoMatrix, Matrix,
                            VertDistance, from.RhoCP, from.RadCP, RadS,
                            SpeciesName, SubCompartName, to.SubCompartName, ScaleName, Test, Regional_and_Continental_deepocean){
  
  # Sedimentation at global scales goes from sea to deeopocean to marinesediment
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & SubCompartName == "sea" & to.SubCompartName == "marinesediment") {
    return(NA)
  }
  
  # If Regional_and_Continental_deepocean is TRUE, the sedimentation rate from sea to deepocean and deepocean to marinesediment
  # at Regional and Continental scale should be NA, because the sedimentation rate goes directly from sea to marinesediment
  # at these scales.  
  if ((ScaleName %in% c("Regional", "Continental")) &&
      ((SubCompartName == "deepocean" && to.SubCompartName == "marinesediment") || (SubCompartName == "sea" && to.SubCompartName == "deepocean")) &&
      (isFALSE(Regional_and_Continental_deepocean) || is.na(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "FALSE")) {
    return(NA)
  }
  
  # If Regional_and_Continental_deepocean is TRUE, remove the sedimentation rate from sea to marinesediment at Regional and Continental scale
  # because deepocean is added between sea and deeopocean.
  if ((ScaleName %in% c("Regional", "Continental")) && SubCompartName == "sea" && to.SubCompartName == "marinesediment" &&
      (Regional_and_Continental_deepocean == "TRUE" || isTRUE(Regional_and_Continental_deepocean))) {
    return(NA)
  }
  
  # No lake and river present at global scales
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & SubCompartName %in% c("lake","river")) {
    return(NA)
  }
  
  switch(SpeciesName,
         "Molecular" = {
           if (to.SubCompartName == "deepocean") {
             return(NA)
           } 
           if (as.character(Test) == "TRUE") {
             if (to.SubCompartName == "lakesediment") {
               return(NA)
             } else {
               SetlingVelocityCP <- 2.5/(24*3600)
               return(SetlingVelocityCP*(1 - FRinw) / VertDistance)
             }
           }
           SetlingVelocityCP <- 
             f_SetVelWater(Shortest_side=from.RadCP*2, 
                           rho_species=from.RhoCP, 
                           rhoMatrix=rhoMatrix, 
                           DynViscWaterStandard=DynViscWaterStandard,
                           DynViscAirStandard=NA,
                           Matrix=Matrix,SubCompartName=SubCompartName, 
                           Shape=NA,
                           Longest_side=NA, Intermediate_side=NA,
                           DragMethod="Original")
           
           return(SetlingVelocityCP*(1 - FRinw) / VertDistance)
         },
         { 
           if (SettlingVelocity <= 0 | is.na(SettlingVelocity)) { # The is.na part was added because for some particles when MinSettVel is na the setvel is not recalculated for some reason which results in an NA value and thus an error in the if statement. 
             return(0)
           }
           return(SettlingVelocity/VertDistance)
         }
         
  )
  
}