#' @title RainOnFreshwater
#' @name RainOnFreshwater
#' @description The fraction of rain falling on lake or river water.
#' @param RAINrate m.s-1
#' @param Area in m2
#' @param SubCompartName #only for lake/rivers
#' @return waterflow of rain directly on lake/river
#' and continental being a part of Moderate (/ Tropic)
#' @export
RainOnFreshwater <- function (RAINrate, Area, SubCompartName, parent, SpeciesName) {
  
  out <- Area |>
    filter(SubCompart %in% c("river", "lake")) |>
    full_join(RAINrate, by="Scale") |>
    expand_grid(SpeciesName) |>
    parent$states$clipStates() |>
    mutate(
      RainOnFreshwater = RAINrate * Area
    ) |>
    filter(!is.na(RainOnFreshwater)) |>
    arrange(Scale, SubCompart) |>
    select(Scale, SubCompart, RainOnFreshwater)
  
  return(data.frame(out))
  
  # if (SubCompartName %in% c("river", "lake")) {
  #   # RAINrateToSI is generarted from units !
  #   return(RAINrate * Area)
  # } else    return(NA)
}
