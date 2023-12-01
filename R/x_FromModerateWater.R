#' @title Mixing of upper and deeper sea water layers (w2 and w3)
#' @name x_FromModerateWater
#' @param OceanCurrent [m3.s-1]
#' @param SubCompartName name of the subcompartment of the box at hand
#' @param ScaleName name of the scale of the box at hand
#' @return Advection sea - deepocean
#' @export
#' 
x_FromModerateWater <- function (to.Volume, to.TAUsea, 
                                    all.x_RegSea2Cont, OceanCurrent, SubCompartName, to.ScaleName) {
  switch(to.ScaleName,
         "Tropic" = {
           switch (SubCompartName,
                   "sea" = {return(0)},
                   "deepocean" = {OceanCurrent}, #fix, data in code
                   NA
           )},
         "Arctic" = {
           switch (SubCompartName,
                   "sea" = {OceanCurrent}, #fix, data in code
                   "deepocean" = {return(0)},
                   NA
           )},
         "Continental" = {
           switch (SubCompartName,
                   "sea" = {
                     RegSea2Cont <- all.x_RegSea2Cont$flow[all.x_RegSea2Cont$fromSubCompart == "sea" & 
                                                             all.x_RegSea2Cont$fromScale == "Regional"]
                     return((to.Volume/to.TAUsea)-RegSea2Cont)}, 
                   NA
           )},
         return(NA)
  )
  
}
