#' @title Mixing of upper and deeper sea water layers (w2 and w3)
#' @name k_AdvectionSeaOcean
#' @param to.Volume Volume of the to box
#' @param from.Volume Volume of the from box
#' @param to.TAUsea TAU - Residence time in sea; scale variable
#' @param from.TAUsea TAU - Residence time in sea; for the "to" box 
#' @param OceanCurrent [m3.s-1]
#' @param SubCompartName name of the subcompartment of the box at hand
#' @param ScaleName name of the scale of the box at hand
#' @return Advection sea - deepocean
#' @export
k_AdvectionSeaOcean <- function (to.Volume, from.Volume,
                     to.TAUsea, from.TAUsea, #either one is NA the other not
                     OceanCurrent, SubCompartName, ScaleName) {

    if (SubCompartName == "sea") {
      
      OceanMixingFlow <- from.Volume / from.TAUsea
      
      if (ScaleName %in% c("Moderate", "Arctic")){ 
        return((OceanMixingFlow + OceanCurrent) / from.Volume)
      } else if (ScaleName == "Tropic") {
        return (OceanMixingFlow / from.Volume)
      } 
      
    } else if (SubCompartName == "deepocean") { #defensive programming; the only other case
      
      OceanMixingFlow <- to.Volume / to.TAUsea
      
      if (ScaleName %in% c("Moderate", "Tropic")){ 
        return((OceanMixingFlow + OceanCurrent) / to.Volume)
      } else if (ScaleName == "Arctic") {
        return (OceanMixingFlow / to.Volume)
      } 
    } 
    
  #all other cases 
  return(NA)
  
}