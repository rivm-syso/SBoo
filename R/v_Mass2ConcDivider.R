#' @title calculate concentration from mass
#' @name Mass2Conc
#' @description calculates Concentrations from mass per Box [g/kg] or [kg/m3]
#' @param Volume of the Box (state) [m3]
#' @param rhoMatrix density of the matrix [kg m-3]
#' @param Fracs fraction of solids in the compartment [-]
#' @param Fracw fraction of water in compartment [-]
#' @return Concentration, kg/kg wet for soils, kg/kg wet for sediment, otherwise kg/m3 
#' @export
Mass2ConcDivider <- function (Volume, Matrix, all.rhoMatrix, Fracs, Fracw) {
  #Some SubComparts do not exist?; have Volume NA
  if (is.na(Volume)) return (NA)
  
  f_Soil.wetweight <- function(Conc.soil, # in kg/m3 soil or sediment
                               Fracw, Fracs){
    Conc.soil*1000/(Fracw*RHOw+Fracs*RHOsolid) # in g/kg (wet) soil
  }

  if (Matrix %in% c("air", "water")) {
    return(1/Volume)
  } else {
    RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "othersoil"]
    RHOw <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "lake"]
    return(f_Soil.wetweight(1 / Volume, Fracw, Fracs, RHOsolid))
  }
}

