#' @title Aerosol water - air Partition coefficient
#' @name Kaerw
#' @description pm
#' @param Kacompw  air water partitioning coefficient [-]
#' @param FRorig fraction of original species [-]
#' @param SubCompartName subcompartment considered
#' @return Kaerw
#' @export
Kaerw <- function (Kacompw, FRorig, SubCompartName, parent, ScaleName, SpeciesName) {
  
  out <- ScaleName |>
    expand_grid(SubCompartName, SpeciesName) |>
    full_join(FRorig, by="SubCompart") |>
    full_join(Kacompw, by="Scale") |>
    parent$states$clipStates() |>
    mutate(
      Kaerw = dplyr::case_when(
        SubCompartName == 'air' ~ 1/(Kacompw*FRorig),
        TRUE ~ NA_real_
      )
    ) |>
    filter(!is.na(Kaerw)) |>
    arrange(Scale, SubCompart) |>
    select(Scale, SubCompart, Kaerw)
  
  return(data.frame(out))
   
  # switch(SubCompartName,
  #        "air" = 1/(Kacompw*FRorig),
  #        NA)
}