#' @title FRACTION of chemical in gas phase of air
#' @name FRingas
#' @description FRACTION of chemical in gas phase of air at EQUILIBRIUM
#' @param FRACw Fraction aerosol water in air
#' @param FRACs Fraction aerosol solids in air
#' @param Kaerw Dimensionless aerosol water / air Partitioning Coefficient
#' @param Kaers Dimensionless aerosol solid / air Partitioning Coefficient
#' @return The fraction of a chemical in the aerosol gas phase. Total: FRingas + FRinaerw + FRinaers = 1.
#' @seealso [Fringas(), FRinw(), FRins()]
#' @export
FRingas <- function(FRACw, FRACs, Kaerw, Kaers, ScaleName, SubCompartName, parent, ...){ #, FRcldw
  
  out <- ScaleName |>
    tidyr::expand_grid(SubCompartName) |>
    parent$states$clipStates(NoSpeciesKey=TRUE) |>
    dplyr::left_join(FRACw, by = c("Scale", "SubCompart")) |>
    dplyr::left_join(FRACs, by = c("Scale", "SubCompart")) |>
    dplyr::left_join(Kaerw, by = c("Scale", "SubCompart")) |>
    dplyr::left_join(Kaers, by = c("SubCompart")) |>
    dplyr::mutate(
      FRingas = 1-FRACw*Kaerw/(1+FRACw*Kaerw+FRACs*Kaers) -FRACs*Kaers/(1+FRACw*Kaerw+FRACs*Kaers)
    ) |>
    dplyr::filter(!is.na(FRingas)) |>
    dplyr::arrange(Scale, SubCompart, Species) |>
    dplyr::select(Scale, SubCompart, Species, FRingas)
  
  return(data.frame(out))
  
  # 1-FRACw*Kaerw/(1+FRACw*Kaerw+FRACs*Kaers) -FRACs*Kaers/(1+FRACw*Kaerw+FRACs*Kaers)
}

