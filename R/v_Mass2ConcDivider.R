#' @title calculate concentration from mass
#' @name Mass2Conc
#' @description calculates Concentrations from mass per Box [g/kg] or [kg/m3]
#' @param Volume of the Box (state) [m3]
#' @param rhoMatrix density of the matrix [kg m-3]
#' @param Fracs fraction of solids in the compartment [-]
#' @param Fracw fraction of water in compartment [-]
#' @return Concentration, kg/kg wet for soils, kg/kg wet for sediment, otherwise kg/m3 
#' @export
Mass2ConcDivider <- function (Volume, Matrix, all.rhoMatrix, FRACa, FRACw) {
  #browser()
  RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "othersoil"]
  RHOw <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "lake"]
  
  ifelse (Matrix %in% c("air", "water"),
    1/Volume,
    { 
      Fracs <- (1-FRACw-FRACa)
      conc <- 1/Volume
      conc/(Fracs * RHOsolid)
    }
  )
}
