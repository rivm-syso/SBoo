#' @title Calculate Area of each compartment
#' @name Area
#' @param AreaLand the land area in the considered compartment [m2]
#' @param AreaSea  the sea area in the considered compartment [m2]
#' @param landFRAC fraction of land in the considered compartment [-]
#' @param SubCompartName the subcompartment which is considered
#' @param ScaleName the scale which is considered
#' @description 
#' Area for air and cloadwater is set to total area of Land and Sea.
#' Area of sea, marinesediment and deepocean is set to area of Sea
#' This is due to the TotalArea being divided in Sea and Land.
#' Area of Land based compartments are all based on the landFRAC (also freshwater compartments).
#' landFRAC and FracSea are variables that are commonly changed in combination with TotalArea to change to a new landscape scenario.
#' Take care of the advective flows as well as they often also need to changed based on the new landscape scenario.
#' 
#' @return Area [m2]
#' @export
Area <- function (AreaLand,
                  AreaSea,
                  landFRAC,
                  all.landFRAC, 
                  SubCompartName,
                  ScaleName) {
  # easiest
  if (SubCompartName %in% c("air", "cloudwater")) {
    return(AreaLand + AreaSea)
  }
  
  if (SubCompartName %in% c("sea", "marinesediment")) {
    return(AreaSea)
  }
  if (SubCompartName == "deepocean" &
      ScaleName %in% c("Arctic", "Moderate", "Tropic")) {
    return(AreaSea)
  }
  
  # if (SubCompartName == "lakesediment" & ScaleName %in% c("Regional", "Continental")){
  #   return(all.landFRAC$landFRAC[all.landFRAC$SubCompart == "lake" & all.landFRAC$Scale == ScaleName] *AreaLand)
  # }
  # 
  # if (SubCompartName == "freshwatersediment" & ScaleName %in% c("Regional", "Continental")){
  #   return(all.landFRAC$landFRAC[all.landFRAC$SubCompart == "river" & all.landFRAC$Scale == ScaleName] *AreaLand)
  # }
  # 
  
  # on land, lake, freshwater and soils;
  # if (ScaleName %in% c("Regional", "Continental")) {
  return(landFRAC * AreaLand)
  # }
  
  # on land, on global scales: only naturalsoil
  # if (SubCompartName == "naturalsoil") {
  #   return(AreaLand)
  # }
  
  #all other cases
  # return(NA) --> These are NA as landFRAC data is NA for all others.
}
