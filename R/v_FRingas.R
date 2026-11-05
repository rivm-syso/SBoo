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
FRingas <- function(FRACw, FRACs, Kaerw, Kaers, parent){ #, FRcldw
  
  out <- FRACw |>
    dplyr::left_join(FRACs) |>
    parent$states$clipStates(NoSpeciesKey=TRUE) |>
    dplyr::left_join(Kaerw) |>
    dplyr::left_join(Kaers) |>
    dplyr::mutate(
      FRingas = 1-FRACw*Kaerw/(1+FRACw*Kaerw+FRACs*Kaers) -FRACs*Kaers/(1+FRACw*Kaerw+FRACs*Kaers)
    )
  
  keepcolumns <- colnames(out)[colnames(out) %in% c(The3D, "FRingas")]
  return(out |> dplyr::select(all_of(keepcolumns)))
  
  # 1-FRACw*Kaerw/(1+FRACw*Kaerw+FRACs*Kaers) -FRACs*Kaers/(1+FRACw*Kaerw+FRACs*Kaers)
}

