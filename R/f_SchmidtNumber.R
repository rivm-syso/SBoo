#' @title SchmidtNumber 
#' @name f_SchmidtNumber
#' @description a description follows .. Joris?
#' @param visc 
#' @param rhoMatrix 
#' @param Temp 
#' @param rad_particle 
#' @return SchmidtNumber  
#' @export
f_SchmidtNumber <- function(visc, rhoMatrix, Diffusivity = NULL, 
                           Temp = NULL, rad_particle = NULL, Cunningham = NULL) {
  stopifnot(is.numeric(visc), is.numeric(rhoMatrix))
  if (is.null(Diffusivity)) {
    if (is.null(Temp) | is.null(rad_particle)) stop("cannot calculate Diffusivity in fSchmidtNumber()")
    if (is.null(Cunningham)) Cunningham <- f_Cunningham(rad_particle)
    Diffusivity <- f_Diffusivity(Temp, visc, rad_particle, Cunningham)
  }
  visc/(rhoMatrix*Diffusivity)
}
