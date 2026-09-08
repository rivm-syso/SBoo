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
    tidyr::expand_grid(SubCompartName, SpeciesName) |>
    dplyr::full_join(FRorig, by="SubCompart") |>
    dplyr::full_join(Kacompw, by="Scale") |>
    parent$states$clipStates() |>
    dplyr::mutate(
      Kaerw = dplyr::case_when(
        SubCompartName == 'air' ~ 1/(Kacompw*FRorig),
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(Kaerw)) |>
    dplyr::arrange(Scale, SubCompart) |>
    dplyr::select(Scale, SubCompart, Kaerw)
  
  return(data.frame(out))
   
  # switch(SubCompartName,
  #        "air" = 1/(Kacompw*FRorig),
  #        NA)
}