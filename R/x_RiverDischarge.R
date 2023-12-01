#' @title RiverDischarge
#' @name x_RiverDischarge
#' @param RainRunoff see ppt
#' @param dischargeFRAC Fraction discharge between freshwater compartment (w1) at regional and continental scale
#' @inheritParams x_ContRiver2Reg we need all.x_ContRiver2Reg to obtain a single number, but
# this way its still (obviously) a flow, to enter tests for fluxes
#' @return River Discharge 
#' @export
x_RiverDischarge <- function (all.Runoff, RainOnFreshwater, 
                              dischargeFRAC, all.x_ContRiver2Reg, 
                              ScaleName, SubCompartName){

  x_ContRiver2Reg <- sum(all.x_ContRiver2Reg$flux) #sum to force an atomic number ?
  SumRainRunoff <- sum(all.Runoff$Runoff[all.Runoff$Scale == ScaleName])
  
    switch (SubCompartName,
    "river" = {
      if(ScaleName == "Continental"){
        
        return((SumRainRunoff + RainOnFreshwater) * (1-dischargeFRAC))
      } 
      if(ScaleName == "Regional"){
        return((SumRainRunoff + RainOnFreshwater + x_ContRiver2Reg) * (1-dischargeFRAC))
      } 
      else NA
    },
    NA
  )
  

}
