#' @title RainOnFreshwater
#' @name RainOnFreshwater
#' @param RAINrate m.s-1
#' @param Area in m2
#' @param SubCompartName #only for lake/rivers
#' @return waterflow of rain directly on lake/river
#' and continental being a part of Moderate (/ Tropic)
#' @export
RainOnFreshwater <- function (RAINrate, Area, SubCompartName) {
  if (SubCompartName %in% c("river", "lake")) {
    # #TODO resolve lake issues in waterflow; for now: old formulas
    # if (SubCompartName == "lake"){
    #   return(0)
    # } else {
    #   # RAINrateToSI is generarted from units !
    #   return(RAINrate * Area)
    # }     
    return(RAINrate * Area)
  } else    return(NA)
}
