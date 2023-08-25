#' @title Diffusivity
#' @name f_Diffusivity
#' @description Diffusivity
#' @param Matrix as in c("air", "water", "soil")
#' @param Temp temperatures
#' @param DynVisc Dynamic viscosity (air, water)
#' @param rad_species see rad_species variable
#' @param Cunningham, see f_Cunningham
#' @return Diffusivity
#' @export
f_Diffusivity <- function(Matrix, Temp, DynVisc, rad_species, Cunningham = NULL) {
  stopifnot(is.numeric(Temp), is.numeric(DynVisc), is.numeric(rad_species))
  if (Matrix == "air") {
    if (is.null(Cunningham))
      Cunningham <- f_Cunningham(rad_species)
    kboltz <- constants::syms$k
    return (kboltz*Temp*Cunningham)/(6*pi*DynVisc*rad_species)
  } else {
    (kboltz*Temp)/(6*pi*visc*rad_species)
  }

}
