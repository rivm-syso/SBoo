#' @title Mixing of upper and deeper sea water layers (w2 and w3)
#' @name x_TropSea2Moder
#' @param OceanCurrent [m3.s-1]
#' @param SubCompartName name of the subcompartment of the box at hand
#' @param ScaleName name of the scale of the box at hand
#' @return Advection sea - deepocean
#' @export
#' 
x_TropSea2Moder <- function (OceanCurrent, SubCompartName, ScaleName) {
  
  switch (SubCompartName,
          "sea" = {OceanCurrent},
          "deepocean" = {return(0)}, #fix, data in code
          NA
  )
  
}