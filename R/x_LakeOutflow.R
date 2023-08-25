#' @title LakeOutflow
#' @name x_LakeOutflow
#' @description related to what flows from the river into the sea
#' @param RiverDischarge see ppt
#' @param ContRiver2Reg see ppt
#' @param LakeFracRiver see ppt
#' @param ScaleName different formula for scale Continental or Regional
#' @return River Discharge 
#' @export
x_LakeOutflow <- function (all.x_RiverDischarge, all.x_ContRiver2Reg, LakeFracRiver, ScaleName){
  x_RiverDischarge <- all.x_RiverDischarge$flow[all.x_RiverDischarge$fromScale==ScaleName]
  if(ScaleName == "Continental"){
    return(LakeFracRiver * x_RiverDischarge)
  } 
  if(ScaleName == "Regional"){
    x_ContRiver2Reg <- all.x_ContRiver2Reg$flow
    return(LakeFracRiver * x_RiverDischarge + x_ContRiver2Reg)
  } 
  else NA
}
