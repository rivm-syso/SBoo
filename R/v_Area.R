#' @title Area
#' @name Area
#' @param AREAland   AREA of land in scale
#' @param AREAsea    AREA of sea in scale
#' @param landFRAC   subcompartment fraction of AREAland
#' @param SubCompartName for which Area is calculated
#' @param ScaleName for which Area is calculated
#' @return Area (but not for sediment) based on data for the SubCompartment / Scale
#' @export
Area <- function (AreaLand,
                     AreaSea,
                     landFRAC,
                     SubCompartName,
                     ScaleName) {
  # easiest
  if (SubCompartName %in% c("air", "cloudwater")) {
    return(AreaLand + AreaSea)
  }
  
  if (SubCompartName == "sea") {
    return(AreaSea)
  }
  
  if (SubCompartName == "deepocean" &
      ScaleName %in% c("Arctic", "Moderate", "Tropic")) {
    return(AreaSea)
  }
  
  # on land, lake, freshwater and soils;
  if (ScaleName %in% c("Regional", "Continental")) {
    return(landFRAC * AreaLand)
  }
  
  # on land, on global scales: only othersoil
  if (SubCompartName == "othersoil") {
    return(AreaLand)
  }
  
  #all other cases
  return(NA)
}
