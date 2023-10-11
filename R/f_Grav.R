#' @title Gravitational impaction collision frequency coefficient (20181102)
#' @name fGrav
#' @description Collission frequency of ENPs with other particulates due to Gravitational or Intertial impaction in s-1 for heteroagglomeration
#' @param viscosity Dynamic viscosity of liquid  (fraction of) compartment [kg.m-1.s-1]
#' @param from.radius Radius of nanoparticle [m]
#' @param from.rho Density of nanoparticle [kg.m-3]
#' @param radius_Otherparticle  Radius of Other particle [m]
#' @param rho_Otherparticle Density (specific weight) of natural particle [kg/m3]
#' @return fGrav [s-1]
#' @export
fGrav <- function(viscosity,from.radius,from.rho,radius_Otherparticle ,rho_Otherparticle,rho.fluid){
  SetVel <- fSetvel(rad_particle=from.radius, #radius particle
                    rho_particle=from.rho, #density particle
                    rho.fluid=rho.fluid,
                    viscosity) #viscosity air
  
  SetVelOther <- fSetvel(rad_particle=radius_Otherparticle , #radius particle
                         rho_particle=rho_Otherparticle, #density particle
                         rho.fluid=rho.fluid,
                         viscosity) #viscosity air
  
  pi*(from.radius+radius_Otherparticle )^2*ABS(SetVel-SetVelOther)
}
