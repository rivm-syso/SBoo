#' @title General Advection process
#' @name k_Advection
#' @description Calculation of k, given a Flow
#' @param flow advection rate [m/s]
#' @param Volume volume of compartment [m3]
#' @param Remove_global If this variable is TRUE, the global scales (Arctic, Moderate and Tropic) are removed
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return Rate constant for 1st order process associated with fluxes
#' @export
#'

k_Advection <- function(flow,
                        Volume,
                        ScaleName,
                        SubCompartName,
                        to.SubCompartName,
                        SpeciesName,
                        Remove_global,
                        Regional_and_Continental_deepocean,
                        AdvInput) {
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) &&
      (!is.na(Remove_global) &&
       (isTRUE(Remove_global) || Remove_global == "TRUE"))) {
    return(NA)
  }
  if(SpeciesName == "Macro"){return(NA)}
  
  # Ensure no flows from regional sea and deepocen to continental deepocean, and no continental river to regional river.
  else if ((
    #(ScaleName == "Regional" && (SubCompartName %in% c("sea", "deepocean")) && to.SubCompartName == "deepocean") ||
    (
      ScaleName == "Regional" &&
      (SubCompartName %in% c("sea")) &&
      to.SubCompartName == "deepocean"
    ) || #there could transfer from deepocean to deepocean
    (
      ScaleName == "Continental" &&
      SubCompartName == "river" && to.SubCompartName == "river"
    )
  ) &&
  (
    !is.na(Regional_and_Continental_deepocean) &&
    (
      isTRUE(Regional_and_Continental_deepocean) ||
      Regional_and_Continental_deepocean == "TRUE"
    )
  )) {
    return(0)
  }
  
  else {
    #if a certain value is input to replace the default, then use it even if it =0
    if (!is.null(AdvInput) &&
        !is.na(AdvInput) && !is.nan(AdvInput)) {
      return(AdvInput)
    }
    
    return(flow / Volume) #not compartment "air"; not a valid airflow
  }
}