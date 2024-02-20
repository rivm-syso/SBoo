#' @title calculate concentration from mass
#' @name Mass2Conc
#' @description calculates Concentrations from mass per Box
#' @param eqMass mass, as returned by solvers
#' @param Volume of the Box (state)
#' @param RHOsolid
#' @param Fracs
#' @param Fracw
#' @return Concentration, kg/kg wet for soils, kg/kg wet for sediment, otherwise kg/m3 
#' @export
Mass2Conc <- function (EqMass, Volume, Matrix, all.rhoMatrix, Fracs, Fracw) {
  
  f_Soil.wetweight <- function(Conc.soil, # in kg/m3 soil or sediment
                               Fracw, Fracs,
                               RHOsolid){
    Conc.soil*1000/(Fracw*RHOw+Fracs*RHOsolid) # in g/kg (wet) soil
  }
  
  # f_Soil.dryweight <- function(Conc.soil, # in kg/m3 soil
  #                              Fracs,
  #                              RHOsolid){
  #   Conc.soil*1000/(Fracs*RHOsolid) # in g/kg (dry) soil
  # }
  
  if (Matrix %in% c("air", "water")) {
    return(EqMass / Volume)
  } else {
    RHOsolid <- all.rhoMatrix$rhoMatrix[SubCompart == "othersoil"]
    RHOw <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "lake"]
    return(f_Soil.wetweight(EqMass / Volume, Fracw, Fracs, RHOsolid))
  }
}

