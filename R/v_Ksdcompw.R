#'@title Dimensionless sed/water PARTITION COEFFICIENT for molecular species specific to sediment compartment
#'@name Ksdcompw
#'@description Dimensionless sed/water PARTITION COEFFICIENT
#'@param FRACw
#'@param FRACs
#'@param Kp (sed/water PARTITION COEFFICIENT) of sediment
#'@param RHOsolid
#'@return 
#'@export
Ksdcompw <- function(FRACw, FRACs, Kp, all.rhoMatrix, Matrix){
  if (Matrix == "sediment") {
    RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
    return(FRACw+FRACs*Kp*RHOsolid/1000) # we need to take care that RHOsolid here can be specific to the compartment compared to generic one used in Kp in relation to Ksw!
  } else
    return(NA)
}
