#'@title FRACa
#'@name FRACa 
#'@description Calculates the total air fraction in the considered system.
#'either subFRACa or, when the main matrix, remainder after substracting subFRACs + subFRACw
#'@param subFRACa subfraction of air in a non-air compartment [-]
#'@param subFRACw subfraction of water in a non-water compartment [-]
#'@param subFRACs subfraction of solids in a non-soil, non-sediment compartment [-]
#'@param Matrix type of compartment 
#'@return FRACa 
#'@export
FRACa <- function(subFRACa, subFRACw, subFRACs, Matrix, parent, ScaleName, SpeciesName) {
  
  out <- ScaleName |>
    expand_grid(SpeciesName, Matrix) |>
    full_join(subFRACw, by=c("Scale", "SubCompart")) |>
    full_join(subFRACs, by=c("Scale", "SubCompart")) |>
    full_join(subFRACa, by=c("Scale", "SubCompart")) |>
    parent$states$clipStates() |>
    mutate(
      FRACa = dplyr::case_when(
        Matrix == 'air' ~ 1 - subFRACw - subFRACs,
        Matrix == 'sediment' ~ 0,
        TRUE ~ subFRACa
      )
    ) |>
    filter(!is.na(FRACa)) |>
    arrange(Scale, SubCompart) |>
    select(Scale, SubCompart, FRACa)
  
  return(data.frame(out))  
  
  # if (Matrix == "air") {
  #   return (1 - subFRACw - subFRACs)
  # } else if (Matrix == "sediment") {
  #   return (0)
  # } else {
  #   return (subFRACa)
  # }
}


