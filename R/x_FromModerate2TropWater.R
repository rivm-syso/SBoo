#' @title Mixing of upper and deeper sea water layers (w2 and w3)
#' @name x_FromModerate2TropWater
#' @param OceanCurrent [m3.s-1]
#' @param SubCompartName name of the subcompartment of the box at hand
#' @param ScaleName name of the scale of the box at hand
#' @return Advection sea - deepocean
#' @export
#' 
x_FromModerate2TropWater <- function ( OceanCurrent, SubCompartName, ScaleName) {
           switch (SubCompartName,
                   "sea" = {return(0)},
                   "deepocean" = {OceanCurrent}, #fix, data in code
                   NA
           )
  
}
