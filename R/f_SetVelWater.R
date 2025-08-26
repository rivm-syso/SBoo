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
f_SetVelWater <- function(rad_species, rho_species, rhoMatrix, 
                     DynViscWaterStandard,
                     DynViscAirStandard,
                     Matrix,SubCompartName, Shape,
                     Longest_side, Intermediate_side, Shortest_side, DragMethod) {
  if (anyNA(c(rho_species,rhoMatrix))){
    return(NA)
  }
  if (is.na(Longest_side) || is.null(Longest_side) || is.na(Intermediate_side) || is.null(Intermediate_side) || is.na(Shortest_side) || is.null(Shortest_side)) {
    Longest_side <- rad_species * 2
    Intermediate_side <- rad_species * 2
    Shortest_side <- rad_species * 2
  } #TODO check application of this default assumption
  if (is.na(Shape) || is.null(Shape)){
    Shape <- "Default"
  }
  
  GN <- constants::syms$gn
  
  if(Matrix == "soil" | Matrix == "sediment") return(NA)
  if(SubCompartName == "cloudwater") return(NA)
  
  if(DragMethod == "Original" & Matrix =="water"){
    sv <- 2*(rad_species^2*(rho_species-rhoMatrix)*GN) / (9*DynViscWaterStandard)
    if (sv <= 0){
      return(0)
    } 
    else {
      return(sv)
    }
  } 
  if(DragMethod == "Original" & Matrix =="air") {
    Cunningham <- f_Cunningham(rad_species)
    sv <- 2*(rad_species^2*(rho_species-rhoMatrix)*GN*Cunningham) / (9*DynViscAirStandard)
    if (sv <= 0){
      return(0)
    } 
    else {
      return(sv)
    }
  } 
  volume <- fVol(rad_species, Shape, Longest_side, Intermediate_side, Shortest_side)
  d_eq <- ( 6/ pi * volume)^(1/3)
  surfaceareaparticle <- f_SurfaceArea(Shape, Longest_side, Intermediate_side, Shortest_side, rad_species)
  surfaceareaperfectsphere <- f_SurfaceArea("Sphere", d_eq, d_eq, d_eq, rad_species)
  #circularity <- Longest_side*Intermediate_side / (d_eq*d_eq)
  perimeterparticle <- f_PerimeterParticle(Shape, Longest_side, Intermediate_side, Shortest_side, rad_species)
  perimetercircle <- f_PerimeterParticle("Sphere", d_eq, d_eq, d_eq, rad_species)
  circularity <- perimeterparticle/perimetercircle
  sphericity <- surfaceareaperfectsphere/surfaceareaparticle
  Psi <- sphericity/circularity # Shape factor Dioguardi
  CSF <- Shortest_side/(sqrt(Longest_side*Intermediate_side)) #Corey Shape Factor
  switch (Matrix,
          "water" = { 
            v_s <- f_SetVelSolver(d_eq=d_eq, Psi=Psi, 
                                  DynViscFluidStandard=DynViscWaterStandard, 
                                  rhoParticle=rho_species, 
                                  rhoFluid=rhoMatrix, DragMethod=DragMethod, 
                                  CSF=CSF, Matrix=Matrix, rad_species=rad_species)
            return(v_s)
          }, 
          "air"= {
            v_s <- f_SetVelSolver(d_eq=d_eq, Psi=Psi, 
                                  DynViscFluidStandard=DynViscAirStandard, 
                                  rhoParticle=rho_species, 
                                  rhoFluid=rhoMatrix, DragMethod=DragMethod, 
                                  CSF=CSF, Matrix=Matrix, rad_species=rad_species)
            return(v_s)
          },
          NA
  )
  if (v_s <= 0) {
    return(0)
  } 
  
  else {
    return(v_s)
  }
}
