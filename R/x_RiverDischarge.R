#' @title RiverDischarge
#' @name x_RiverDischarge
#' @param RainRunoff see ppt
#' @param dischargeFRAC see ppt
#' @param ContRiver2Reg see ppt
#' @return River Discharge 
#' @export
x_RiverDischarge <- function (all.Runoff, RainOnFreshwater, dischargeFRAC, all.x_ContRiver2Reg, ScaleName){
  # twisted, exceptional; we need all.x_ContRiver2Reg to obtain a single number, but\
  # this way its still (obiously) a flow, to enter tests for fluxes
  # rain runoff from soils and direct rain on river
  x_ContRiver2Reg <- sum(all.x_ContRiver2Reg$flux) #sum to force an atomic number ?
  SumRainRunoff <- sum(all.Runoff$Runoff[all.Runoff$Scale == ScaleName])
  if(ScaleName == "Continental"){
    return((SumRainRunoff + RainOnFreshwater) * (1-dischargeFRAC))
  } 
  if(ScaleName == "Regional"){
    return((SumRainRunoff + RainOnFreshwater + x_ContRiver2Reg) * (1-dischargeFRAC))
  } 
  else NA
}
