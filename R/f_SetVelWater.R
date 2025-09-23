#' @title Function for calculating the Settling Velocity of a particle through water
#' @name f_SetVelWater
#' @description Settling Velocity or terminal Velocity of particles in a fluid medium based on Stokes' theorem (1901)
#' NB SettlingVelocity is a function to be used for all types of partical species, not just nano-species. Therefor it's not a variable defining function.
#' @param rhoParticle Density of particle [kg/m3]
#' @param rhoWater Density of fluid matrix in which particle is present [kg/m3]
#' @param DynViscWaterStandard Dynamic viscosity of water []
#' @param radius Radius of the particle [m]
#' @param GN gravitational force constabt [m2/s]
#' @return f_SetVelWater
#' @export
f_SetVelWater <- function(Shortest_side, # for simplification only Shortest_side is needed (is thus 2 x radius!)
                          rho_species, rhoMatrix,
                          DynViscWaterStandard,
                          DynViscAirStandard,
                          Matrix, SubCompartName, Shape,
                          Longest_side, Intermediate_side, DragMethod) {
  if (anyNA(c(rho_species, rhoMatrix))) {
    return(NA)
  }
  if (Matrix == "soil" | Matrix == "sediment") {
    return(NA)
  }
  if (SubCompartName == "cloudwater") {
    return(NA)
  }
  # Check if any of Intermediate or Longest sides is NA or NULL and assign default values if so
  if (is.na(Intermediate_side) || is.null(Intermediate_side) || is.na(Longest_side) || is.null(Longest_side)) {
    Intermediate_side <- Shortest_side * 2 #maybe 0.75 or build in shape functions
    Longest_side <- Shortest_side * 2 # maybe consider other default based on shape
  }
  rad_particle = Shortest_side/2 # for Originial and Cunningham in air settling rate

  GN <- constants::syms$gn
  
  if (is.na(Shape) || is.null(Shape)) {
    Shape <- "Default"
  }

  if (DragMethod == "Original" & Matrix == "water") {
    sv <- 2 * (rad_particle^2 * (rho_species - rhoMatrix) * GN) / (9 * DynViscWaterStandard)
    if (sv <= 0) {
      return(0)
    } else {
      return(sv)
    }
  }
  if (DragMethod == "Original" & Matrix == "air") {
    Cunningham <- f_Cunningham(rad_particle)
    sv <- 2 * (rad_particle^2 * (rho_species - rhoMatrix) * GN * Cunningham) / (9 * DynViscAirStandard)
    if (sv <= 0) {
      return(0)
    } else {
      return(sv)
    }
  }
  
  volume <- fVol(Shape=Shape, Longest_side=Longest_side, Intermediate_side=Intermediate_side, Shortest_side=Shortest_side)
  d_eq <- (6 / pi * volume)^(1 / 3)
  surfaceareaparticle <- f_SurfaceAreaParticle(Shape=Shape, Longest_side=Longest_side, Intermediate_side=Intermediate_side, Shortest_side=Shortest_side)
  surfaceareaperfectsphere <- f_SurfaceAreaParticle(Shape="Sphere", rad_species=d_eq/2)
  # circularity <- Longest_side*Intermediate_side / (d_eq*d_eq)
  perimeterparticle <- f_PerimeterParticle(Shape=Shape, Longest_side=Longest_side, Intermediate_side=Intermediate_side, Shortest_side=Shortest_side)
  perimetercircle <- f_PerimeterParticle(Shape="Sphere", rad_species=d_eq/2)
  circularity <- perimeterparticle / perimetercircle
  sphericity <- surfaceareaperfectsphere / surfaceareaparticle
  Psi <- sphericity / circularity # Shape factor Dioguardi
  CSF <- Shortest_side / (sqrt(Longest_side * Intermediate_side)) # Corey Shape Factor
  #Parameters for Bagheri et al. 2016
  alpha <- 0.45+10/exp(2.5*log10(rho_species/rhoMatrix)+30) 
  beta <- 1-37/exp(3*log10(rho_species/rhoMatrix)+100)
  f <- Shortest_side/Intermediate_side
  e <-  Intermediate_side/Longest_side
  FN <- f^2*e*(d_eq^3/(Longest_side*Intermediate_side*Shortest_side))
  FS <- f*e^1.3*(d_eq^3/(Longest_side*Intermediate_side*Shortest_side))
  kS <- 1/2*(FS^(1/3)+FS^(-1/3))
  kN <- 10^(alpha*(-log10(FN))^beta)
  
  switch(Matrix,
    "water" = {
      v_s <- f_SetVelSolver(
        d_eq = d_eq, Psi = Psi,
        DynViscFluidStandard = DynViscWaterStandard,
        rhoParticle = rho_species,
        rhoFluid = rhoMatrix, DragMethod = DragMethod,
        CSF = CSF, Matrix = Matrix, rad_species = NA, kS, kN
      )
      return(v_s)
    },
    "air" = {
      v_s <- f_SetVelSolver(
        d_eq = d_eq, Psi = Psi,
        DynViscFluidStandard = DynViscAirStandard,
        rhoParticle = rho_species,
        rhoFluid = rhoMatrix, DragMethod = DragMethod,
        CSF = CSF, Matrix = Matrix, rad_species = rad_particle, kS, kN
      )
      return(v_s)
    },
    NA
  )
  if (v_s <= 0) {
    return(0)
  } else {
    return(v_s)
  }
}
