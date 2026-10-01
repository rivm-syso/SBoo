#'@title Fraction of intercepted particles in soil from Runoff
#'@name IntrcptFrac 
#'@description Data from the interception experiments in Runoff from Han et al. 2022 are used to 
#'
#'@param subFRACa fraction of air in a non-air subcompartment [-]
#'@param subFRACw fraction of water in a non-water subcompartment [-]
#'@param subFRACs fraction of solids in a non-soil, non-sediment compartment [-]
#'@param Matrix type of compartment
#'@return FRACs 
#'@export
IntrcptFrac <- function(VegInterceptFrac, SizeRunoff, rad_species, Matrix, parent, SpeciesName){
  
  out <- rad_species |>
    dplyr::left_join(VegInterceptFrac, by=c("SubCompart")) |>
    dplyr::left_join(SizeRunoff, by=c("Scale")) |>
    dplyr::filter(SubCompart %in% c("othersoil", "naturalsoil", "agriculturalsoil")) |>
    parent$states$clipStates() |>
    dplyr::mutate(
      interceptFraction = dplyr::case_when(
        rad_species < SizeRunoff ~ 0,
        TRUE ~ VegInterceptFrac
      )
    ) |>
    dplyr::filter(!is.na(interceptFraction)) |>
    dplyr::select(Scale, SubCompart, Species, interceptFraction) |>
    dplyr::arrange(Scale, SubCompart, Species)
  
  return(data.frame(out))
    
  # if (Matrix %in% c("soil")) {
  #   interceptFraction = ifelse(rad_species < SizeRunoff,  
  #                              0, 
  #                              VegInterceptFrac)
  #   return(interceptFraction)
  # }
}

