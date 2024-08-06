#' @title calculate concentration from mass
#' @name Mass2Conc
#' @description calculates Concentrations from mass per Box [g/kg] or [kg/m3]
#' @param eqMass mass, as returned by solvers [kg]
#' @param Volume of the Box (state) [m3]
#' @param rhoMatrix density of the matrix [kg m-3]
#' @param Fracs fraction of solids in the compartment [-]
#' @param Fracw fraction of water in compartment [-]
#' @return Concentration, kg/kg wet for soils, kg/kg wet for sediment, otherwise kg/m3 
#' @export
Mass2Conc <- function (EqMass, Volume, Matrix, all.rhoMatrix, Fracs, Fracw) {
  
  f_Soil.wetweight <- function(Conc.soil, # in kg/m3 soil or sediment
                               Fracw, Fracs,
                               RHOsolid){
    Conc.soil*1000/(Fracw*RHOw+Fracs*RHOsolid) # in g/kg (wet) soil
  }

  if (Matrix %in% c("air", "water")) {
    return(EqMass / Volume)
  } else {
    RHOsolid <- all.rhoMatrix$rhoMatrix[SubCompart == "othersoil"]
    RHOw <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "lake"]
    return(f_Soil.wetweight(EqMass / Volume, Fracw, Fracs, RHOsolid))
  }
}

