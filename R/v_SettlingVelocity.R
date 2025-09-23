#' @title SettlingVelocity
#' @name v_SettlingVelocity
#' @description Settling Velocity or terminal Velocity of particles in a fluid medium based on different DragMethods.
#' NB SettlingVelocity is a function to be used for all types of partical species, not just nano-species. Therefor it's not a variable defining function.
#' @param rho_species Density of particle [kg/m3]
#' @param rhoMatrix Density of fluid matrix in which particle is present [kg/m3]
#' @param DynViscAirStandard Dynamic viscosity of the fluid matrix [kg m-1 s-1]
#' @param DynViscWaterStandard Dynamic viscosity of the fluid matrix [kg m-1 s-1]
#' @param rad_species Radius of the particle [m]
#' @param Matrix function is defined for Water and Air; slightly different algorithm [text]
#' @param SubCompartName name of the different subcompartments, used to segregate different formulas [text] 
#' @param Shape Shape as defined by user, different possibilities (see f_Vol for options) [text]
#' @param Longest_side the longest side of the particle as defined by the user [m]
#' @param Intermediate_side the intermediate side of the particle as defined by user [m]
#' @param Shortest_side the shortst side of the particle as identified by user [m]
#' @param DragMethod The Method used for computing the drag coefficient as defined by user, opportunity for choosing 4 different ones. See f_DragCoefficient for options
#' @return Settling velocity [m.s-1]
#' @export
SettlingVelocity <- function(rad_species, rho_species, rhoMatrix, 
                             DynViscWaterStandard,
                             DynViscAirStandard,
                             Matrix,SubCompartName, ScaleName,
                             Shape,Longest_side,
                             Intermediate_side, DragMethod) {
  if (anyNA(c(rho_species,rhoMatrix))){
    return(NA)
  }
  if ((ScaleName %in% c("Regional", "Continental")) & SubCompartName == "deepocean") {
    return(NA)
  }
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & SubCompartName %in% c("lake","river")) {
    return(NA)
  }
  # Check if Shortest side is NA or NULL and assign default values if so
  # if ( is.na(Shortest_side) || is.null(Shortest_side) ) {
  #   Shortest_side <- rad_particle * 2
  # }
  # Check if any of Intermediate or Longest sides is NA or NULL and assign default values if so

  
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
  
  if (is.na(Intermediate_side) || is.null(Intermediate_side) ||is.na(Longest_side) || is.null(Longest_side)) {
    Intermediate_side <- rad_species * 2 #maybe 0.75 or build in shape functions
    Longest_side <- rad_species * 2
    warning("Need for Intermediate_side or Longest_side, but not provided, setting to 2*rad_species")
  }
  
  Species_Volume <- fVol(rad_particle = rad_species,
                         Shape = Shape, 
                         Longest_side = Longest_side,
                         Intermediate_side = Intermediate_side)
  
  d_eq <- (6/ pi * Species_Volume)^(1/3) # calculate equivalent diameter of perfect sphere
  
  surfaceareaparticle <- f_SurfaceAreaParticle(Shape=Shape, 
                                               Intermediate_side=Intermediate_side, 
                                               Longest_side=Longest_side, 
                                               rad_particle=rad_species)
  surfaceareaperfectsphere <- f_SurfaceAreaParticle(Shape="Sphere", rad_particle=d_eq/2)
  #circularity <- Longest_side*Intermediate_side / (d_eq*d_eq)
  perimeterparticle <- f_PerimeterParticle(Shape=Shape, Intermediate_side=Intermediate_side,Longest_side=Longest_side, rad_particle=rad_species)
  perimetercircle <- f_PerimeterParticle(Shape="Sphere", rad_particle=d_eq/2)
  circularity <- perimeterparticle/perimetercircle
  sphericity <- surfaceareaperfectsphere/surfaceareaparticle
  Psi <- sphericity/circularity # Shape factor Dioguardi
  CSF <- rad_species/(sqrt(Longest_side*Intermediate_side)) #Corey Shape Factor
  #Parameters for Bagheri et al. 2016
  alpha <- 0.45+10/exp(2.5*log10(rho_species/rhoMatrix)+30) 
  beta <- 1-37/exp(3*log10(rho_species/rhoMatrix)+100)
  f <- Shortest_side/Intermediate_side
  e <-  Intermediate_side/Longest_side
  FN <- f^2*e*(d_eq^3/(Longest_side*Intermediate_side*Shortest_side))
  FS <- f*e^1.3*(d_eq^3/(Longest_side*Intermediate_side*Shortest_side))
  kS <- 1/2*(FS^(1/3)+FS^(-1/3))
  kN <- 10^(alpha*(-log10(FN))^beta)
  
  switch (Matrix,
          "water" = { 
            v_s <- f_SetVelSolver(d_eq=d_eq, Psi=Psi, 
                                  DynViscFluidStandard=DynViscWaterStandard, 
                                  rhoParticle=rho_species, 
                                  rhoFluid=rhoMatrix, DragMethod=DragMethod, 
                                  CSF=CSF, Matrix=Matrix, rad_species=rad_species,
                                  kS=kS, kN=kN)
            return(v_s)
          }, 
          "air"= {
            v_s <- f_SetVelSolver(d_eq=d_eq, Psi=Psi, 
                                  DynViscFluidStandard=DynViscAirStandard, 
                                  rhoParticle=rho_species, 
                                  rhoFluid=rhoMatrix, DragMethod=DragMethod, 
                                  CSF=CSF, Matrix=Matrix, rad_species=rad_species,
                                  kS=kS, kN=kN)
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

