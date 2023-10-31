#' @title Brownian motion or kinetic collision frequency coefficient (20181102)
#' @name f_Brown
#' @description Collission frequency of ENPs with other particulates due to bownian motion (kinetic energy) in s-1 for heteroagglomeration
#' @param Temp temperature [K]
#' @param viscosity Dynamic viscosity of liquid  (fraction of) compartment [kg.m-1.s-1]
#' @param from.radius Radius of nanoparticle [m]
#' @param radius_Otherparticle  Radius of Other particle [m]
#' @return fBrown [s-1]
#' @export
f_Brown <- function(Temp,viscosity,from.radius,radius_Otherparticle){
  ((2*getConst("r")*Temp)/(3*viscosity))*(from.radius+radius_Otherparticle )^2/(from.radius*radius_Otherparticle)
}
