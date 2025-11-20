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
k_Sedimentation <- function(FRinw, SettlingVelocity, DynViscWaterStandard,
                            VertDistance, from.RhoCP, from.RadCP, RadS,
                            SpeciesName, SubCompartName, to.SubCompartName, ScaleName, Test, Test_surface_water){
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & SubCompartName == "sea" & to.SubCompartName == "marinesediment") {
    return(NA)
  }
  
  #if Test_surface_water is TRUE, then compute a ksed for deepocean and sea at regional and continental scale
  if ((ScaleName %in% c("Regional", "Continental")) &&
      (SubCompartName == "deepocean") &&
      (isFALSE(Test_surface_water) || is.na(Test_surface_water))) {
    return(NA)
  }
  
  #if Test_surface_water is TRUE, remove sea -> marinesediment (replaced by sea -> deepocean -> marinsediment)
  if ((ScaleName %in% c("Regional", "Continental")) && SubCompartName == "sea" && to.SubCompartName == "marinesediment" &&
      isTRUE(Test_surface_water)) {
    return(NA)
  }
  if ((ScaleName %in% c("Regional", "Continental")) & SubCompartName == "deepocean") {
    return(NA)
  }
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
           if (SettlingVelocity <= 0) {
             return(0)
           }
           return(SettlingVelocity/VertDistance)
         }
         
  )
  
}