#' @title Mixing of upper and deeper sea water layers (w2 and w3)
#' @name x_OceanMixing2Sea
#' @param to.Volume Volume of the to box
#' @param from.Volume Volume of the from box
#' @param to.TAUsea TAU - Residence time in sea; scale variable
#' @param from.TAUsea TAU - Residence time in sea; for the "to" box 
#' @param OceanCurrent [m3.s-1]
#' @param SubCompartName name of the subcompartment of the box at hand
#' @param ScaleName name of the scale of the box at hand
#' @return Advection sea - deepocean
#' @export
x_OceanMixing2Sea <- function (all.Volume,
                               TAUsea,
                               OceanCurrent, SubCompartName, ScaleName) {
  
  
  
  switch (SubCompartName,
          "deepocean" = {
            toVolume <- all.Volume$Volume[all.Volume$SubCompart == "sea" &
                                            all.Volume$Scale == ScaleName]
           OceanMixingFlow <- toVolume / TAUsea
            
            if (ScaleName %in% c("Moderate", "Tropic", "Arctic")){ 
              return((OceanMixingFlow + OceanCurrent) )
            } else return (NA)
          },
          return(NA)
  )
  
}

