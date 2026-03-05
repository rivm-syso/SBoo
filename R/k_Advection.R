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

k_Advection <- function(flow, Volume, ScaleName, SubCompartName, to.SubCompartName,to.ScaleName, Remove_global, Regional_and_Continental_deepocean, AdvInput) { 
  if (ScaleName %in% c("Tropic", "Moderate", "Arctic") & (!is.na(Remove_global) && (isTRUE(Remove_global) || Remove_global == "TRUE"))) {
    return(NA) }
  
  # Ensure no flows from regional sea and deepocen to continental deepocean, and no continental river to regional river. 
  else if (
    ((ScaleName == "Regional" & (SubCompartName == "sea" | SubCompartName == "deepocean")) |
      ((ScaleName == "Continental" & SubCompartName == "river")) & (to.ScaleName == "Regional" & to.SubCompart == "river")) &
    (!is.na(Regional_and_Continental_deepocean) && (isTRUE(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "TRUE"))
  ) {
    return(0)
    
  } else { 

    #if a certain value is input to replace the default, then use it
    if (!is.null(AdvInput) && !is.na(AdvInput) && !is.nan(AdvInput)) {
      return(AdvInput)
    }
  
    return(flow/Volume) #not compartment "air"; not a valid airflow 
  } 
}
