#' @title Fraction of chemical in aerosol water
#' @name FRinaerw
#' @description Calculates the fraction of a chemical in aerosol water relative to the total amount of the chemical in air.
#' @param Kaers Dimensionless aerosol solid / air Partitioning Coefficient [-]
#' @param Kaerw Dimensionless aerosol water / air Partitioning Coefficient [-]
#' @param FRACs fraction of solids in the matrix [-]
#' @param FRACw fraction of water in the matrix [-]
#' @returns The fraction of a chemical in the aerosol water phase. Total: FRingas + FRinaerw + FRinaers = 1.
#' @seealso [Fringas(), FRinw(), FRins()]
#' @export
FRinaerw <- function (Kaers, Kaerw, FRACw, FRACs, parent) {

  out <- FRACw |>
    dplyr::full_join(FRACs, by = c("Scale", "SubCompart")) |>
    dplyr::full_join(Kaerw, by = c("Scale", "SubCompart")) |>
    dplyr::full_join(Kaers, by = "SubCompart") |>
    parent$states$clipStates(NoSpeciesKey=TRUE) |>
    dplyr::mutate(
      FRinaerw = FRACw*Kaerw/(1+FRACw*Kaerw+FRACs*Kaers)
    ) |>
    dplyr::filter(!is.na(FRinaerw)) |>
    dplyr::arrange(Scale, SubCompart) |>
    dplyr::select(Scale, SubCompart, FRinaerw)
    
  return(data.frame(out))
  
  # FRACw*Kaerw/(1+FRACw*Kaerw+FRACs*Kaers) 
}

