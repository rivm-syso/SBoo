#Script for calculating Run-off and Erosion

#'@title Run-off of particles (unbound to be added)
#'@name k_Runoff
#'@param FRACrun Fraction run-off #[-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance Mixing depth soil #[m]
#'@param EROSIONsoil Soil erosion #[mm/yr]
#'@return k_Run-off Run-off of particles from soil to water #[s-1]
#'@export

k_Runoff <- function(FRACrun, RAINrate, VertDistance, Kscompw,
                     relevant_depth_s, penetration_depth_s,
                     ScaleName, to.SubCompartName, SpeciesName){
  if (ScaleName %in% c("Regional", "Continental") & to.SubCompartName == "sea") {
    return(NA)
  } 
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName != "sea") {
    return(NA)
  } 
  switch(SpeciesName,
         "Molecular" = {
           ( RAINrate * FRACrun / Kscompw)* 
             f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s) / VertDistance #[s-1]
         },
         (FRACrun*RAINrate*
            f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s))/VertDistance
         #Runoff for particulates should maybe be made dependent on the size of particles that flows freely with water? As opposed to erosion? A function of shear, with a bariere basedon certrain high intensity rain episode?
  )
}
