#' @title SchmidtNumber 
#' @name f_SchmidtNumber
#' @description Calculation of the dimensionless Schmidt number, based on Wang et al. (2010) https://doi.org/10.5194/acp-10-5685-2010
#' Denotes the ratio of momentum diffusivity and mass diffusivity.
#' @param visc Dynamic viscosity of air [kg.m-1.s-1]
#' @param rhoMatrix Density of the compartment the particle is in [kg m-3]
#' @param Temp Temperature of the compartment the particle is in [kg m-3]
#' @param rad_particle radius of the particle [m]
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
