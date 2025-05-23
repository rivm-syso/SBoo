#'@title Run-off of particles
#'@name k_Runoff
#'@description Runoff rate corrected for soil penetration depth for molecules and particulates
#'@param Runoff Runoff, as computed by v_Runoff [s-1]
#'@param FracROWatComp fraction of water component runoff can flow to [-]
#'@param Volume volume of compartment 
#'@param Kscompw soil water partitioning coefficient [-]
#'@param relevant_depth_s the soil depth at which a process occurs [m]
#'@param penetration_depth_s assumed penetration depth from Hollander et al. (2007), see f_CORRsoil [m]
#'@param ScaleName considered ScaleName
#'@param SubCompartName subcompartment considered
#'@param SpeciesName considered species 
#'@param Matrix type of supcompartment 
#'@param FRACrun Fraction of precipation of run-off #[-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance Mixing depth soil #[m]
#'@return k_Run-off Run-off of particles from soil to water #[s-1]
#'@export

k_Runoff <- function(Runoff,VertDistance, to.FracROWatComp,
                     Volume, Kscompw,
                     relevant_depth_s, penetration_depth_s,
                     ScaleName, to.SubCompartName, to.ScaleName, SpeciesName, Matrix, all.Matrix){
  if (ScaleName %in% c("Regional", "Continental") & to.SubCompartName == "sea") {
    return(NA)
  } 
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName != "sea") {
    return(NA)
  } 

  switch(SpeciesName,
         "Molecular" = {
           (Runoff / Kscompw) * 
             f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s) / Volume  * to.FracROWatComp  #[s-1]
         },
         (Runoff *
            f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s))/ Volume * to.FracROWatComp 
         
  )
}
