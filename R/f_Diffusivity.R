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
  if(!is.numeric(Temp)){
    warning("Temp missing in f_Diffusivity")
    return(NA)
  }
  if(!is.numeric(DynVisc)){
    warning("DynVisc missing in f_Diffusivity")
    return(NA)
  }
  if(!is.numeric(rad_species)){
    warning("rad_species missing in f_Diffusivity")
    return(NA)
  }
  kboltz <- constants::syms$k
  if (Matrix == "air") {
    if (is.null(Cunningham))
      Cunningham <- f_Cunningham(rad_species)
    return ((kboltz*Temp*Cunningham)/(6*pi*DynVisc*rad_species))
  } else {
    (kboltz*Temp)/(6*pi*DynVisc*rad_species)
  }

}
