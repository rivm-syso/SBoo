#' @title Fraction of chemical in aerosol water
#' @name FRinaerw
#' @description Calculates the fraction of a chemical in aerosol water relative to the total amount of the chemical in air.
#' @param FRACaer Fraction aerosols in air [-]
#' @inheritParams Kaers
#' @inheritParams Kaerw
#' @returns The fraction of a chemical in the aerosol water phase. Total: FRingas + FRinaerw + FRinaers = 1.
#' @seealso [Fringas(), FRinw(), FRins()]
#' @export
FRinaerw <- function (Kaers, Kaerw, FRACaerw, FRACaers) {
  FRACaerw*Kaerw/(1+FRACaerw*Kaerw+FRACaers*Kaers) 
}

