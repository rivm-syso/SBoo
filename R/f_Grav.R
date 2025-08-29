#' @title Gravitational impaction collision frequency coefficient (20181102)
#' @name f_Grav
#' @description Collission frequency of ENPs with other particulates due to Gravitational or Intertial impaction in s-1 for heteroagglomeration
#' @param DynViscWaterStandard Dynamic viscosity of liquid  (fraction of) compartment [kg.m-1.s-1]
#' @param Shortest_side Diameter (2 x Radius) or Shortest side of the particle (nanomaterial or microplastics) [m]
#' @param Intermediate_side Intermediate side of the particle (nanomaterial or microplastics) [m]
#' @param Longest_side Longest side of the particle (nanomaterial or microplastics) [m]
#' @param rho Density of nanoparticle [kg.m-3]
#' @param radius_Otherparticle  Radius of Other particle [m]
#' @param rho_Otherparticle Density (specific weight) of natural particle [kg/m3]
#' @return f_Grav [s-1]
#' @export
f_Grav <- function(rhoParticle,
                   radius_Otherparticle, rho_Otherparticle,
                   rhoFluid, DynViscWaterStandard,
                   Matrix, SubCompartName, Shape,
                   Shortest_side,
                   Longest_side, Intermediate_side, DragMethod) {
  
  SetVel <- f_SetVelWater(
    rho_species = rhoParticle,
    rhoMatrix = rhoFluid,
    DynViscWaterStandard = DynViscWaterStandard,
    DynViscAirStandard = NA,
    Matrix = Matrix, SubCompartName = SubCompartName,
    Shape = Shape,
    Longest_side = Longest_side,
    Intermediate_side = Intermediate_side,
    Shortest_side = Shortest_side,
    DragMethod = DragMethod
  )

  SetVelOther <- f_SetVelWater(
    rho_species = rho_Otherparticle,
    rhoMatrix = rhoFluid,
    DynViscWaterStandard,
    DynViscAirStandard = NA,
    Matrix = Matrix, SubCompartName = SubCompartName,
    Shape = NA,
    Longest_side = NA,
    Intermediate_side = NA,
    Shortest_side = Shortest_side,
    DragMethod = "Original"
  )


  pi * (Shortest_side/2 + radius_Otherparticle)^2 * abs(SetVel - SetVelOther)
}
