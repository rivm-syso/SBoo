#' @title Aerosol water - air Partition coefficient
#' @name Kaerw
#' @description Fraction of some fraction
#' @param Kaw 
#' @param Corg
#' @param Rho
#' @param FRorig
#' @param FRACaer
#' @return FRears[]
#' @export
Kaerw <- function (Kacompw, FRorig, Matrix) {

  switch(Matrix,
         "air" = 1/(Kacompw*FRorig),
         NA)
  

  
}