#' @title Fraction of chemical in aerosol solids
#' @name FRinaers
#' @description Calculates the fraction of a chemical in aerosol solids relative to the total amount of the chemical in air.
#' @param FRACaer Fraction aerosols in air [-]
#' @inheritParams Kaers
#' @inheritParams Kaerw
#' @returns The fraction of a chemical in the aerosol solid phase. Total: FRingas + FRinaerw + FRinaers = 1.
#' @seealso [Fringas(), FRinw(), FRins()]
#' @export
FRinaers <- function (Kaers, Kaerw, FRACaers, FRACaerw) {
  FRACaers*Kaers/(1+FRACaerw*Kaerw+FRACaers*Kaers)
}
