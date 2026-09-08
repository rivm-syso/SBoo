#'@title FRACw
#'@name FRACw
#'@description Fraction of water in any matrix, 
#' either subFRACw or, when the main matrix, remainder after substracting subFRACs + subFRACa
#'@param subFRACa subfraction of air in a non-air compartment [-]
#'@param subFRACw subfraction of water in a non-water compartment [-]
#'@param subFRACs subfraction of solids in a non-soil, non-sediment compartment [-]
#'@param Matrix type of compartment 
#'@return FRACw
#'@export
FRACw <- function(subFRACa, subFRACw, subFRACs, Matrix, parent, ScaleName, SpeciesName){
  
  out <- ScaleName |>
    tidyr::expand_grid(SpeciesName, Matrix) |>
    dplyr::full_join(subFRACw, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(subFRACa, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(subFRACs, by=c("Scale", "SubCompart")) |>
    parent$states$clipStates() |>
    dplyr::mutate(
      FRACw = dplyr::case_when(
        Matrix == 'water' ~ 1 - subFRACs - subFRACa,
        TRUE ~ subFRACw
      )
    ) |>
    dplyr::filter(!is.na(FRACw)) |>
    dplyr::arrange(Scale, SubCompart) |>
    dplyr::select(Scale, SubCompart, FRACw)

  return(data.frame(out))  
  
  # if (Matrix == "water") {
  #   return (1 - subFRACs - subFRACa)
  # } else
  #   return (subFRACw)
}
