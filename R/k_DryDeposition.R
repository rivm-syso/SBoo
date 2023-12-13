#' @title Dry deposition from air to surface soil or water
#' @name k_DryDeposition
#' @description 
#' Calculation of the first order rate constant for deposition from air to the soil or water surface [s-1]
#' Diffusivity is calculated from f_Diffusivity
#' @param to.Area Surface area of receiving land or water compartment [m2]
#' @param from.Volume Volume of air compartment [m3]
#' @param Temp Temperature [K]
#' @param rhoMatrix Density of air [kg.m-3]
#' #param from.viscosity Dynamic viscosity of air compartment [kg.m-1.s-1]
#' @param ColRad Land surface particle collector radius [m] appendix A of LOTEUR.
#' #param FricVel Friction velocity [m s-1] (19 according to Van Jaarsveld, table 2.2)
#' @param to.alpha.surf depends on the vegetation type, see table A.1 in ref guide LOTEUR v2.0 2016
#' @param AEROresist aerodynamic resistance (Constant set to 74) [s.m-1]
#' @param DynViscAirStandard Dynamic viscosity of air
#' @param rad_species rad_particle
#' @param rho_species rho_particle
#' @param gamma.surf
#' @param from.SettlingVelocity Settlings Velocity of particulate species in air [m.s-1]
#' @param Matrix
#' @return k_drydeposition, the rate constant for 1rst order process: dry deposition from air to soil or water [s-1]
#' @export
k_DryDeposition <- function(to.Area, from.Volume, AEROresist, to.gamma.surf, FricVel,
                            DynViscAirStandard, rhoMatrix, to.ColRad, rad_species, rho_species,
                            Temp,to.alpha.surf, from.SettlingVelocity, SpeciesName, SubCompartName){
  
  if (SpeciesName %in% c("Nanoparticle","Aggregated","Attached")) {
    if (anyNA(c(AEROresist, DynViscAirStandard, rhoMatrix, rho_species, to.alpha.surf))) {
      return(NA)
    }
    switch(SubCompartName,
           "air" = {
             
             Cunningham <- f_Cunningham(rad_species)
             Diffusivity <- f_Diffusivity(Matrix = "air", Temp, DynViscAirStandard, rad_species, Cunningham)
             
             SchmidtNumber <- DynViscAirStandard/(rhoMatrix*Diffusivity) #rhoMatrix to be converted to RhoWater or RhoAir
             # alpha.surf = depends on vegetation type, see e.g. LOTEUR ref guide table A.1
             Brown <- SchmidtNumber^(-to.gamma.surf)
             
             StN <- ifelse(to.ColRad==0|is.na(to.ColRad), # StokesNumber following ref guide LOTEUR v2.0 2016
                           (from.SettlingVelocity*FricVel)/DynViscAirStandard/rhoMatrix, # for smooth surfaces (water)
                           (from.SettlingVelocity*FricVel)/(constants::syms$gn*to.ColRad)      # for vegetated surfaces (soil)
             )
             Intercept <- ifelse(to.ColRad==0|is.na(to.ColRad),
                                 0, # for smooth surfaces
                                 0.5 * (rad_species/to.ColRad)^2 # LOTEUR (eq 5.14) in ref guide.
             )
             
             #StokesNumberVeg <- fStokesNumber_rain(rad_particle, rho_species, from.rho, from.visc,RAINrate,FRACtwet,Cunningham.cw,g) # to be update in future!
             R1 <- exp(-StN^0.5) # R1 = correction factor representing the fraction of particles that stick to the surface
             epsilon = 3 # epsilon  is empirical constant set to 3
             beta.a = 2 # constant set to 2
             
             # in SB4N alpha.surf is 0.8!!!
             Impaction <- (StN/(StN+to.alpha.surf))^beta.a # LOTUS EUROS ref guide (eq. 5.13)  
             
             SurfResist <- 1/(epsilon*FricVel*(Brown+ # Collection efficiency for Brownian diffusion
                                                 Intercept+ # Collection efficiency for interception
                                                 Impaction)*R1) # Collection efficiency for impaction
             DRYDEPvelocity <- 1/(AEROresist+SurfResist)+from.SettlingVelocity
             
             # Currently implemented in SimpleBox for P species:
             # if (species == "P"){
             #   DRYDEPvelocity <- AEROSOLdeprate # AEROSOLdeprate constant given in xls version of SB4
             # }
             (DRYDEPvelocity*to.Area)/from.Volume
             
           },
           NA)
    
  } else {
    return (NA)
  }
}


