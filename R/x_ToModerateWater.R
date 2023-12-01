#' @title Mixing of upper and deeper sea water layers (w2 and w3)
#' @name x_ToModerateWater
#' @param OceanCurrent [m3.s-1]
#' @param SubCompartName name of the subcompartment of the box at hand
#' @param ScaleName name of the scale of the box at hand
#' @return Advection sea - deepocean
#' @export
#' 
x_ToModerateWater <- function (Volume, TAUsea, 
                             all.x_RegSea2Cont,OceanCurrent, SubCompartName, ScaleName) {
  switch(ScaleName,
         "Tropic" = {
           switch (SubCompartName,
                   "sea" = {OceanCurrent},
                   "deepocean" = {return(0)}, #fix, data in code
                   NA
           )},
         "Arctic" = {
           switch (SubCompartName,
                   "sea" = {return(0)}, #fix, data in code
                   "deepocean" = {OceanCurrent},
                   NA
           )},
         "Continental" = {
           switch (SubCompartName,
                   "sea" = {
                     RegSea2Cont <- all.x_RegSea2Cont$flow[all.x_RegSea2Cont$fromSubCompart == "sea" & 
                                                             all.x_RegSea2Cont$fromScale == "Regional"]
                     return((Volume/TAUsea)-RegSea2Cont)}, 
                   NA
           )},
         return(NA)
  )
  
}
