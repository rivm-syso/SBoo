#' @title Fraction of chemical in aerosol solids
#' @name FRinaers
#' @description Calculates the fraction of a chemical in aerosol solids relative to the total amount of the chemical in air.
#' @param Kaers Dimensionless aerosol solid / air Partitioning Coefficient [-]
#' @param Kaerw Dimensionless aerosol water / air Partitioning Coefficient [-]
#' @param FRACs fraction of solids in the matrix [-]
#' @param FRACw fraction of water in the matrix [-]
#' @returns The fraction of a chemical in the aerosol solid phase. Total: FRingas + FRinaerw + FRinaers = 1.
#' @seealso [Fringas(), FRinw(), FRins()]
#' @export
FRinaers <- function (Kaers, Kaerw, FRACs, FRACw, parent) {
  
  out <- FRACw |>
    dplyr::full_join(FRACs, by = c("Scale", "SubCompart")) |>
    dplyr::full_join(Kaerw, by = c("Scale", "SubCompart")) |>
    dplyr::full_join(Kaers, by = "SubCompart") |>
    parent$states$clipStates(NoSpeciesKey=TRUE) |>
    dplyr::mutate(
      FRinaers = FRACs*Kaers/(1+FRACw*Kaerw+FRACs*Kaers)
    ) |>
    dplyr::filter(!is.na(FRinaers)) |>
    dplyr::arrange(Scale, SubCompart) |>
    dplyr::select(Scale, SubCompart, FRinaers)
  
  return(data.frame(out))
  
  # FRACs*Kaers/(1+FRACw*Kaerw+FRACs*Kaers)
}
