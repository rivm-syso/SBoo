#' @title FRACTION of chemical in gas phase of air
#' @name FRingas
#' @description FRACTION of chemical in gas phase of air at EQUILIBRIUM
#' @param FRACw Fraction aerosol water in air
#' @param FRACs Fraction aerosol solids in ai
#' @param Kaerw Dimensionless aerosol water / air Partitioning Coefficient
#' @param Kaers Dimensionless aerosol solid / air Partitioning Coefficient
#' @return FRgas[]
#' @export
FRingas <- function(FRACw, FRACs, Kaerw, Kaers, ...){ #, FRcldw
  1-FRACw*Kaerw/(1+FRACw*Kaerw+FRACs*Kaers) -FRACs*Kaers/(1+FRACw*Kaerw+FRACs*Kaers)
}

