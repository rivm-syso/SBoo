#' @title AirFlow
#' @name AirFlow
#' @param Volume in m3
#' @param Area in m2
#' @param WINDspeed m.s-1
#' @param SubCompartName for which the calculation is executed
#' @return AirFlow, from one scale to another 
#' @export
AirFlow <- function (Volume, Area, WINDspeed, SubCompartName){
  
  if (SubCompartName %in%  c("air")) { #, "cloudwater" should also?!
    TAU <- f_TAU(Area, WINDspeed) #Residence time
    Volume / TAU
  } else {
    return(NA) #not compartment "air"; not a valid airflow
  }
}
