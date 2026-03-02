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
  if (ScaleName %in% c("Tropic", "Moderate", "Arctic") & (!is.na(Remove_global) && (isTRUE(Remove_global) || Remove_global == "TRUE"))) {
    return(NA) }
  
  # This was done in ScaleSubcompartData before, but why would we want to remove the advection flows from sea and deepocean at Regional scale to Continental scale?
  else if (ScaleName == "Regional" & (!is.na(Test_surface_water) && (isTRUE(Test_surface_water) || Test_surface_water == "TRUE")) ){
    return(0)
    
  } else { 

    #if a certain value is input to replace the default, then use it
    if (!is.null(AdvInput) && !is.na(AdvInput) && !is.nan(AdvInput)) {
      return(AdvInput)
    }
  
    return(flow/Volume) #not compartment "air"; not a valid airflow 
    } }
