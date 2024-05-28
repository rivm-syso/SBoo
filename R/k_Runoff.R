#Script for calculating Run-off and Erosion

#'@title Run-off of particles (unbound to be added)
#'@name k_Runoff
#'@param FRACrun Fraction run-off #[-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance Mixing depth soil #[m]
#'@return k_Run-off Run-off of particles from soil to water #[s-1]
#'@export

k_Runoff <- function(Runoff,VertDistance, FracROWatComp,
                     Volume, Kscompw,
                     relevant_depth_s, penetration_depth_s,
                     ScaleName, to.SubCompartName, to.ScaleName, SpeciesName, Matrix, all.landFRAC, all.Matrix){
  if (ScaleName %in% c("Regional", "Continental") & to.SubCompartName == "sea") {
    return(NA)
  } 
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName != "sea") {
    return(NA)
  } 
  #fraction <- FracROWatComp(all.landFRAC, all.Matrix, Matrix, SubCompartName = to.SubCompartName, ScaleName = ScaleName)
  switch(SpeciesName,
         "Molecular" = {
           (Runoff / Kscompw) * 
             f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s) / Volume  * FracROWatComp  #[s-1]
         },
         (Runoff *
            f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s))/ Volume * FracROWatComp 
         #Runoff for particulates should maybe be made dependent on the size of particles that flows freely with water? As opposed to erosion? A function of shear, with a bariere basedon certrain high intensity rain episode?
  )
}
