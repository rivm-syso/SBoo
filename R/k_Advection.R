#' @title General Advection process
#' @name k_Advection
#' @description Calculation of k, given a Flow
#' @param flow advection rate [m/s]
#' @param Volume volume of compartment [m3]
#' @return Rate constant for 1st order process associated with fluxes
#' @export
#' 
#' 
#' 

k_Advection <- function(flow, Volume, ScaleName, SubCompartName, to.SubCompartName, Remove_global, Test_surface_water, AdvInput) { 
  if (ScaleName %in% c("Tropic", "Moderate", "Arctic") 
      & (!is.na(Remove_global) && (isTRUE(Remove_global) || Remove_global == "TRUE"))) {
    return(NA) } 
  else { 
    #if (ScaleName %in% c("Regional") 
    #    & (!is.na(Test_surface_water) && (isTRUE(Test_surface_water) || Test_surface_water == "TRUE"))) { #Use the k_advection provided by default - Hajjar et al. 2025
   # }
    
    #if a certain value is input to replace the default, then use it
    if (!is.null(AdvInput) && !is.na(AdvInput) && !is.nan(AdvInput) && AdvInput != 0) {
      return(AdvInput)
    }
    
    
    
    return(flow/Volume) #not compartment "air"; not a valid airflow 
    } }
#  if ((fromScale %in% c("Tropic", "Moderate", "Artic")) ||
 #     (toScale   %in% c("Tropic", "Moderate", "Artic")) &&
  #    (!is.na(Remove_global) && 
   #    (isTRUE(Remove_global) || Remove_global == "TRUE"))) { 
    