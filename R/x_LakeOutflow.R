#' @title LakeOutflow
#' @name x_LakeOutflow
#' @description related to what flows from the river into the sea
#' @param RiverDischarge see ppt
#' @param ContRiver2Reg see ppt
#' @param LakeFracRiver Fraction (default 0.1) of r
#' @param ScaleName different formula for scale Continental or Regional
#' @return River Discharge 
#' @export
x_LakeOutflow <- function (all.x_RiverDischarge, 
                           all.x_ContRiver2Reg, 
                           LakeFracRiver, 
                           ScaleName,
                           all.landFRAC){
  x_RiverDischarge <- all.x_RiverDischarge$flow[all.x_RiverDischarge$fromScale==ScaleName]
  if(ScaleName == "Continental"){
    landFRACContinental <- all.landFRAC[all.landFRAC$Scale == "Continental", ]
    LakeAreaFrac <- landFRACContinental$landFRAC[landFRACContinental$SubCompart == "lake"]
    RiverAreaFrac <- landFRACContinental$landFRAC[landFRACContinental$SubCompart == "river"]
    LakeFracRiver <- LakeAreaFrac/(RiverAreaFrac+LakeAreaFrac) #Calculate the % of the river discharge that comes from the lake
    return(LakeFracRiver * x_RiverDischarge)
  } 
  if(ScaleName == "Regional"){
    landFRACContinental <- all.landFRAC[all.landFRAC$Scale == "Regional", ]
    LakeAreaFrac <- landFRACContinental$landFRAC[landFRACContinental$SubCompart == "lake"]
    RiverAreaFrac <- landFRACContinental$landFRAC[landFRACContinental$SubCompart == "river"]
    LakeFracRiver <- LakeAreaFrac/(RiverAreaFrac+LakeAreaFrac)
    x_ContRiver2Reg <- all.x_ContRiver2Reg$flow
    return(LakeFracRiver * (x_RiverDischarge + x_ContRiver2Reg))
  } 
  else NA
}
