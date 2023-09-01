#'@title Dimensionless soil/water Partition Coefficient of molecular species specific to soil subcompartment
#'@name Kscompw
#'@description Dimensionless soil/water PARTITION COEFFICIENT
#'@param FRACw
#'@param FRACa NB FRACs is calculated as the remiander from FRACa and FRACw
#'@param Kp (soil/water PARTITION COEFFICIENT)
#'@param RHOsolid
#'@return 
#'@export
Kscompw <- function(FRACw, FRACa, Kacompw, FRorig_spw, Kp, all.rhoMatrix, Matrix){
  if (Matrix == "soil") {
    RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
    return(FRACa*(Kacompw*FRorig_spw)+FRACw+(1-FRACa-FRACw)*Kp*RHOsolid/1000) 
    # we need to take care that RHOsolid above can be specific to the compartment
    # compared to generic one used in Kp in relation to Ksw!
  } else
    return(NA)
}
