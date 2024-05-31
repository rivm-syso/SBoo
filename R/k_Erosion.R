#'@title Erosion of particles
#'@name k_Erosion
#'@param CORRrunoff correction for runoff
#'@param VertDistance Mixing depth soil #[m]
#'@param EROSIONsoil Soil erosion #[mm/yr]
#'@param ScaleName To adjust for the absence of rivers in global scales
#'@param to.SubCompartName name of the subcompartment of the destination box of this process
#'@return k_Erosion Erosion of attached enp species (P) from soil to water #[s-1]
#'@export
k_Erosion <- function(relevant_depth_s,penetration_depth_s, EROSIONsoil, VertDistance,
                      to.FracROWatComp,
                      ScaleName, to.SubCompartName,  Matrix, all.landFRAC, all.Matrix ){
  if (ScaleName %in% c("Regional", "Continental") & to.SubCompartName == "sea") {
    return(NA)
  } 
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName != "sea") {
    return(NA)
  } 
 # fraction <- FracROWatComp(all.landFRAC, all.Matrix, Matrix, SubCompartName, ScaleName)
  EROSIONsoil * f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s) / VertDistance * to.FracROWatComp  #[s-1]
}
