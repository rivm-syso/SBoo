#'@title Kp
#'@name Kp
#'@description subcompartment/water PARTITION COEFFICIENT
#'@param FRorig, 
#'@param Ksw, 
#'@param Ksw.alt, 
#'@param RHOsolid
#'@param Corg
#'@param CorgSoilStandard
#'@param ChemClass Chemical class, e.g. neutral or metal
#'@param Matrix the medium, the formula is only applicable to soil and sediment
#'@export
Kp <- function(FRorig, Ksw, Ksw.alt, all.rhoMatrix, Corg, CorgStandard, Matrix, ChemClass){
  if (Matrix %in% c("soil", "sediment","water")) {
    RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
    return(
      switch (ChemClass,
        "Acid" = (FRorig*Ksw + (1-FRorig)*Ksw.alt) * (1000 / RHOsolid) * (Corg / CorgStandard),
        "Base" = (FRorig*Ksw + (1-FRorig)*Ksw.alt) * (1000 / RHOsolid) * (Corg / CorgStandard),
        {(FRorig*Ksw) * (1000 / RHOsolid) * (Corg / CorgStandard)} #Corg
      )
      
    )
  } else return (NA)
}
