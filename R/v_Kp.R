#'@title Kp
#'@name Kp
#'@description subcompartment/water PARTITION COEFFICIENT
#'@param FRorig, 
#'@param KswDorC, 
#'@param KswDorC.alt, 
#'@param RHOsolid
#'@param Corg
#'@param CorgSoilStandard
#'@param ChemClass Chemical class, e.g. neutral or metal
#'@param Matrix the medium, the formula is only applicable to soil and sediment
#'@export
Kp <- function(FRorig, KswDorC, Ksw.alt, all.rhoMatrix, Corg, CorgStandard, Matrix, ChemClass){
  if (Matrix %in% c("soil", "sediment","water")) {
    RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
    return(
      switch (ChemClass,
        "acid" = (FRorig*KswDorC + (1-FRorig)*Ksw.alt) * (1000 / RHOsolid) * (Corg / CorgStandard),
        "base" = (FRorig*KswDorC + (1-FRorig)*Ksw.alt) * (1000 / RHOsolid) * (Corg / CorgStandard),
        {(FRorig*KswDorC) * (1000 / RHOsolid) * (Corg / CorgStandard)} #Corg
      )
      
    )
  } else return (NA)
}
