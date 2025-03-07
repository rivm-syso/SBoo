#' @title LakeOutflow
#' @name x_LakeOutflow
#' @description related to what flows from the river into the sea
#' @param RiverDischarge see ppt
#' @param ContRiver2Reg see ppt
#' @param LakeFracRiver Fraction (default 0.1) of r
#' @param ScaleName different formula for scale Continental or Regional
#' @return River Discharge 
#' @export
x_LakeOutflow <- function (RainOnFreshwater,
                           #all.x_RiverDischarge, 
                           #all.x_ContRiver2Reg,
                           #LakeFracRiver, 
                           all.Runoff,
                           FracROWatComp,
                           ScaleName){
  # browser()
  #x_RiverDischarge <- all.x_RiverDischarge$flow[all.x_RiverDischarge$fromScale==ScaleName]
  SumRainRunoff <- sum(all.Runoff$Runoff[all.Runoff$Scale == ScaleName])
  if(ScaleName == "Continental"){
    return(RainOnFreshwater + FracROWatComp*SumRainRunoff)
    # return(LakeFracRiver * x_RiverDischarge)
  } 
  if(ScaleName == "Regional"){
    # x_ContRiver2Reg <- all.x_ContRiver2Reg$flow
    # return(LakeFracRiver * (x_RiverDischarge + x_ContRiver2Reg)) # x_ContRiver2Reg is faulty!
    return(RainOnFreshwater + FracROWatComp*SumRainRunoff)
  } 
  else NA
}
