#' @title SettlingVelocity
#' @name v_SettlingVelocity
#' @description Settling Velocity or terminal Velocity of particles in a fluid medium based on ...(REF)
#' NB SettlingVelocity is a function to be used for all types of partical species, not just nano-species. Therefor it's not a variable defining function.
#' @param rho_species Density of particle [kg/m3]
#' @param rhoMatrix Density of fluid matrix in which particle is present [kg/m3]
#' @param DynViscAirStandard Dynamic viscosity of the fluid matrix [unit]
#' @param DynViscWaterStandard Dynamic viscosity of the fluid matrix [unit]
#' @param rad_species Radius of the particle [m]
#' @param Matrix function is defined for Water and Air; slightly different algorithm [text]
#' @return Settling velocity
#' @export
SettlingVelocity <- function(rad_species, rho_species, rhoMatrix, 
                             DynViscWaterStandard,
                             DynViscAirStandard,
                             Matrix,SubCompartName) {
  if (anyNA(c(rho_species,rhoMatrix))){
    return(NA)
  }

  GN <- constants::syms$gn
  
  switch (Matrix,
          "water" = {
            if(SubCompartName == "cloudwater") return(NA)
            2*(rad_species^2*(rho_species-rhoMatrix)*GN) / (9*DynViscWaterStandard)
            
          },
          "air" = {
            Cunningham <- f_Cunningham(rad_species)
            2*(rad_species^2 * (rho_species - rhoMatrix)*GN*Cunningham)/(9*DynViscAirStandard) # particle settling velocity
          },
          NA #else
  )
}
