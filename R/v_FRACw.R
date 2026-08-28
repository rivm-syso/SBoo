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
    expand_grid(SpeciesName, Matrix) |>
    full_join(subFRACw, by=c("Scale", "SubCompart")) |>
    full_join(subFRACa, by=c("Scale", "SubCompart")) |>
    full_join(subFRACs, by=c("Scale", "SubCompart")) |>
    parent$states$clipStates() |>
    mutate(
      FRACw = dplyr::case_when(
        Matrix == 'water' ~ 1 - subFRACs - subFRACa,
        TRUE ~ subFRACw
      )
    ) |>
    filter(!is.na(FRACw)) |>
    arrange(Scale, SubCompart) |>
    select(Scale, SubCompart, FRACw)

  return(data.frame(out))  
  
  # if (Matrix == "water") {
  #   return (1 - subFRACs - subFRACa)
  # } else
  #   return (subFRACw)
}
