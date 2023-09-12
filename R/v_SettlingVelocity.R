#' @title SettlingVelocity
#' @name v_SettlingVelocity
#' @description Settling Velocity or terminal Velocity of particles in a fluid medium based on ...(REF)
#' NB SettlingVelocity is a function to be used for all types of partical species, not just nano-species. Therefor it's not a variable defining function.
#' @param rho_species Density of particle [kg/m3]
#' @param matrix.Rho Density of fluid matrix in which particle is present [kg/m3]
#' @param DynVisc Dynamic viscosity of the fluid matrix [unit]
#' @param rad_species Radius of the particle [m]
#' @param Matrix function is defined for Water and Air; slightly different algorithm [text]
#' @return Settling velocity
#' @export
SettlingVelocity <- function(rad_species, rho_species, matrix.Rho, DynVisc,
                             Matrix) {
  stopifnot(is.numeric(rho_species), is.numeric(matrix.Rho), is.numeric(DynVisc))

  GN <- constants::syms$gn
  
  switch (Matrix,
          "water" = {
            #     ifelse(ScaleName == "Regional", # This is not logical, should work for all scales
            2*(rad_species^2*(rho_species-matrix.Rho)*GN) / (9*DynVisc) #,
            #        SettlVelocitywater # what is this? Does not make sense.
            #      )
          },
          "air" = {
            Cunningham <- f_Cunningham(rad_species)
            2*(rad_species^2 * (rho_species - matrix.Rho)*GN*Cunningham)/(9*DynVisc) # particle settling velocity
          },
          NA #else
  )
}
