#'@title FRACs
#'@name FRACs 
#'@description fraction of solid in any matrix 
#' either subFRACs or, when the main matrix, remainder after substracting subFRACa + subFRACw
#'@param subFRACa fraction of air in a non-air subcompartment [-]
#'@param subFRACw fraction of water in a non-water subcompartment [-]
#'@param subFRACs fraction of solids in a non-soil, non-sediment compartment [-]
#'@param Matrix type of compartment
#'@return FRACs 
#'@export
FRACs <- function(subFRACa, subFRACw, subFRACs, Matrix, parent, ScaleName, SpeciesName){
  
  out <- ScaleName |>
    expand_grid(SpeciesName, Matrix) |>
    full_join(subFRACw, by=c("Scale", "SubCompart")) |>
    full_join(subFRACs, by=c("Scale", "SubCompart")) |>
    full_join(subFRACa, by=c("Scale", "SubCompart")) |>
    parent$states$clipStates() |>
    mutate(
      subFRACa = dplyr::case_when(
        Matrix == 'sediment' ~ 0,
        TRUE ~ subFRACa
      ),
      FRACs = dplyr::case_when(
        Matrix %in% c("soil", "sediment") ~ 1 - subFRACw - subFRACa,
        TRUE ~ subFRACs
      )
    ) |>
    arrange(Scale, SubCompart) |>
    filter(!is.na(FRACs)) |>
    select(Scale, SubCompart, FRACs)
    
  return(data.frame(out))  

  # if (Matrix %in% c("soil", "sediment")) {
  #   if (Matrix == "sediment") subFRACa <- 0
  #   return (1 - subFRACw - subFRACa)
  # } else
  #   return (subFRACs)
}

