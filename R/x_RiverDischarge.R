#' @title RiverDischarge [s-1]
#' @name x_RiverDischarge
#' @param Runoff Runoff flow from soil to river [m3.s-1]
#' @param RainOnFreshwater Water flow of rain directly on lake/river [m3.s-1] 
#' @param dischargeFRAC Fraction discharge regional fresh water to continental scale [-]
#' @param x_ContRiver2Reg Flow from continental to regional river water [s-1] 
#' @param ScaleName Name of the scale of the box at hand
#' @param SubCompartName Name of the subcompartment of the box at hand
#' @return River Discharge [s-1]
#' @export
x_RiverDischarge <- function (all.Runoff, all.RainOnFreshwater, 
                              dischargeFRAC, all.x_ContRiver2Reg, 
                              ScaleName, SubCompartName){
  x_ContRiver2Reg <- sum(all.x_ContRiver2Reg$flow) #sum to force an atomic number ?
  SumRainRunoff <- sum(all.Runoff$Runoff[all.Runoff$Scale == ScaleName]) +
    sum(all.RainOnFreshwater$RainOnFreshwater[all.RainOnFreshwater$Scale == ScaleName])
  
  switch (SubCompartName,
          "river" = {
            if(ScaleName == "Continental"){
              return((SumRainRunoff) * (1-dischargeFRAC))
            } 
            if(ScaleName == "Regional"){
              return((SumRainRunoff + x_ContRiver2Reg) * (1-dischargeFRAC))
            } 
            else NA
          },
          NA
  )
  
  
}
