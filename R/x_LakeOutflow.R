#' @title LakeOut
#' @name x_LakeOut
#' @description This calculates the flow from a freshwater compartment within the regional or continental scale to the main freshwater compartment (w1).
#' @param RainOnFreshwater The direct rain on the freshwater compartment (e.g. a lake).
#' @param all.Runoff Runoff from soil to the freshwater compartment
#' @param FracROWatComp The fraction runoff discharging into the seperate w0 freshwater compartment based on surface area fraction.
#' @param ScaleName different formula for scale Continental or Regional
#' @return River Discharge 
#' @export
x_LakeOut <- function (RainOnFreshwater,
                           #all.x_RiverDischarge, 
                           #all.x_ContRiver2Reg,
                           #LakeFracRiver, 
                           all.Runoff,
                           FracROWatComp,
                           ScaleName){
  # browser()
  SumRunoff <- sum(all.Runoff$Runoff[all.Runoff$Scale == ScaleName])
  return(RainOnFreshwater + FracROWatComp*SumRunoff)
  #x_RiverDischarge <- all.x_RiverDischarge$flow[all.x_RiverDischarge$fromScale==ScaleName]
  # SumRunoff <- sum(all.Runoff$Runoff[all.Runoff$Scale == ScaleName])
  # if(ScaleName == "Continental"){
  #   return(RainOnFreshwater + FracROWatComp*SumRunoff)
  #   # return(LakeFracRiver * x_RiverDischarge)
  # } 
  # if(ScaleName == "Regional"){
  #   # x_ContRiver2Reg <- all.x_ContRiver2Reg$flow
  #   # return(LakeFracRiver * (x_RiverDischarge + x_ContRiver2Reg))
  #   return(RainOnFreshwater + FracROWatComp*SumRunoff)
  # } 
  # else NA
}
