#' @title Dry deposition from air to surface soil or water of particulate species
#' @name k_DryDeposition
#' @description Calculation of the first order rate constant for deposition from air to the soil or water surface [s-1]
#' Adjustments for (Test - FALSE, new SB version) are made based on LOTOS-EUROS reference guide v.2.2002 https://www.rivm.nl/lotos-euros
#' @param to.Area Surface area of receiving land or water compartment [m2]
#' @param from.Volume Volume of air compartment [m3]
#' @param Temp Temperature [K]
#' @param rhoMatrix Density of air [kg.m-3]
#' #param from.viscosity Dynamic viscosity of air compartment [kg.m-1.s-1]
#' @param ColRad Land surface particle collector radius [m] appendix A of LOTEUR.
#' @param FricVel Friction velocity [m s-1] (19 according to Van Jaarsveld, table 2.2)
#' @param to.alpha.surf depends on the vegetation type, see table A.1 in ref guide LOTEUR v2.0 2016
#' @param AEROresist aerodynamic resistance (Constant set to 74) [s.m-1]
#' @param DynViscAirStandard Dynamic viscosity of air
#' @param rad_species rad_particle [m]
#' @param rho_species rho_particle [kg m3]
#' @param gamma.surf resistance parameter depending on average surface type as taken from the LOTOS-EUROS reference guide table A.1
#' @param rhoMatrix density of the matrix [kg m3]
#' @param SpeciesName name of the species considered
#' @param SubcompartName name of the Subcompartment 
#' @param AEROSOLdeprate deposition velocity of aerosol particles [m/s]
#' @param Test determines if SB4-Excel approach is taken or enhanced method from R version [boolean]
#' @param from.SettlingVelocity Settlings Velocity of particulate species in air [m.s-1]
#' @param Matrix type of compartment [-]
#' @param Cunningham collision frequency of ENPs with other particles [-] see f_Cunnningham
#' @param Diffusivity see f_Diffusivity
#' @param FRACtwet fraction of wet periods [-]
#' @return k_drydeposition, the rate constant for 1rst order process: dry deposition from air to soil or water [s-1]
#' @export
k_DryDeposition <- function(Area, Volume, AEROresist, gamma.surf, FricVel,
                            DynViscAirStandard, rhoMatrix, ColRad, rad_species, rho_species,
                            Temp, alpha.surf, SettlingVelocity, SpeciesName, SubCompartName,
                            AEROSOLdeprate, Test = FALSE, Matrix, RAINrate, FRACtwet, ScaleName, parent) {
  
  rho_water <- 998
  AREAFRACveg <- 0.01
  LargeVegRadius <- 5e-4
  VegHair <- 1e-5
  epsilon <- 3 # epsilon is empirical constant set to 3
  beta.a <- 2 # constant set to 2
  sea_AERO <- 135

  # Names of dataframe
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_DryDeposition")))

  # Build grids, from
  args <- c("ScaleName", "SubCompartName", "SpeciesName", "Matrix", "Volume", "rhoMatrix", "rad_species", "rho_species", "Temp", "SettlingVelocity", "AEROSOLdeprate", "RAINrate", "FRACtwet")
  dfs <- mget(args, envir = environment())
  grid_from <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()

  # To
  args <- c("ScaleName", "SubCompartName", "SpeciesName", "Matrix", "Area", "gamma.surf", "ColRad", "alpha.surf")
  dfs <- mget(args, envir = environment())
  grid_to <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()

  # Output
  out <- parent$FromDataAndTo("k_DryDeposition") |>
    # From
    dplyr::left_join(grid_from, by=c("from.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species")) |>

    #To
    dplyr::left_join(grid_to, by=c("to.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species"), suffix=c("_from", "_to")) |>

    dplyr::mutate(
      Cunningham_rad_species = f_Cunningham(rad_species),
      rad_RainDrop = f_RadRain(RAINrate, FRACtwet),
      Cunningham.cw = f_Cunningham(rad_RainDrop),
      Diffusivity = f_Diffusivity(Matrix = c("air"), Temp, DynViscAirStandard, rad_species, Cunningham_rad_species), #air hard coded, per regel al niet al goed
      SchmidtNumber = DynViscAirStandard / (rhoMatrix * Diffusivity),
      Brown = dplyr::case_when(
        from.SubCompart == 'air' & Matrix_to == 'water' & as.character(Test) == 'TRUE' ~ SchmidtNumber^-(1/2),
        from.SubCompart == 'air' & as.character(Test) == 'TRUE' ~ SchmidtNumber^(-2/3),
        from.SubCompart == 'air' & as.character(Test) == 'FALSE' ~ SchmidtNumber^(-gamma.surf),
        TRUE ~ NA
      ),
      Settvel.Particle.cw = ((2*rad_RainDrop)^2 * (rho_water - rhoMatrix)* #rhoMatrix to be converted to RhoWater and RhoAir
                                 constants::syms$gn*Cunningham.cw)/(18*DynViscAirStandard),
      Relax.Particle.a = ((rho_species-rhoMatrix)*(2*rad_species)^2*Cunningham_rad_species)/(18*DynViscAirStandard),
      StN = dplyr::case_when(
        from.SubCompart == 'air' & as.character(Test) == 'TRUE' ~ (2*Relax.Particle.a*(Settvel.Particle.cw-SettlingVelocity))/(2*rad_RainDrop),
        # StokesNumber following ref guide LOTEUR v2.0 2016
        from.SubCompart == 'air' & as.character(Test) == 'FALSE' & (ColRad == 0 | is.na(ColRad)) ~ (SettlingVelocity * FricVel) / DynViscAirStandard / rhoMatrix, # for smooth surfaces (water)
        from.SubCompart == 'air' & as.character(Test) == 'FALSE' & (ColRad != 0 | !is.na(ColRad)) ~ (SettlingVelocity * FricVel) / (constants::syms$gn * ColRad), # for vegetated surfaces (soil)
        TRUE ~ NA
      ),
      Intercept = dplyr::case_when(
        from.SubCompart == 'air' & as.character(Test) == 'TRUE' ~ 0.3*(AREAFRACveg * (rad_species/(rad_species+VegHair))+(1-AREAFRACveg)*(rad_species/(rad_species+LargeVegRadius))),
        from.SubCompart == 'air' & as.character(Test) == 'FALSE' & (ColRad == 0 | is.na(ColRad)) ~ 0, # for smooth surfaces
        from.SubCompart == 'air' & as.character(Test) == 'FALSE' & (ColRad != 0 | !is.na(ColRad)) ~ 0.5 * (rad_species / ColRad)^2, # LOTEUR (eq 5.14) in ref guide.
        TRUE ~ NA
      ),
      R1 = exp(-StN^0.5),
      # in SB4N alpha.surf is 0.8!!!
      Impaction = dplyr::case_when(
        from.SubCompart == 'air' & as.character(Test) == 'TRUE' & Matrix_to == 'water' ~  10^(-3/StN),
        from.SubCompart == 'air' & as.character(Test) == 'TRUE' & Matrix_to == 'soil' ~ (StN / (StN + 0.8))^beta.a,
        from.SubCompart == 'air' & as.character(Test) == 'FALSE' ~ (StN / (StN + alpha.surf))^beta.a,
        TRUE ~ NA
      ),
      SurfResist = dplyr::case_when(
        from.SubCompart == 'air' & as.character(Test) == 'TRUE' & Matrix_to == 'water' ~  1 /(FricVel * (Brown + Impaction)),
        from.SubCompart == 'air' & as.character(Test) == 'TRUE' & Matrix_to != 'water' ~ 1 / (FricVel * (Brown + # Collection efficiency for Brownian diffusion
                                                                                                           Intercept + # Collection efficiency for interception
                                                                                                           Impaction)),
        from.SubCompart == 'air' & as.character(Test) == 'FALSE' ~ 1 / (epsilon * FricVel * (Brown + # Collection efficiency for Brownian diffusion
                                                                                               Intercept + # Collection efficiency for interception
                                                                                               Impaction) * R1), # Collection efficiency for impaction,
        TRUE ~ NA
      ),
      DRYDEPvelocity = dplyr::case_when(
        from.SubCompart == 'air' & SpeciesName_from != "Attached" ~ 1 / (AEROresist + SurfResist) + SettlingVelocity,
        from.SubCompart == 'air' & SpeciesName_from == "Attached" ~ AEROSOLdeprate,
        TRUE ~ NA
      ),
      k_DryDeposition = dplyr::case_when(
        SpeciesName_from %in% c("Nanoparticle", "Aggregated", "Attached") & from.SubCompart == 'air' ~ (DRYDEPvelocity * Area) / Volume,
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(k_DryDeposition)) |>
    dplyr::select(names, k_DryDeposition) |>
    dplyr::arrange(dplyr::across(names))

  return(data.frame(out))

  
  # 
  # switch(as.character(Test),
  #        "TRUE" = {
  #          # Existing function logic when test is FALSE
  #          if (SpeciesName %in% c("Nanoparticle", "Aggregated", "Attached")) {
  #            if (anyNA(c(AEROresist, DynViscAirStandard, rhoMatrix, rho_species, to.alpha.surf))) {
  #              return(NA)
  #            }
  #            switch(SubCompartName,
  #                   "air" = {
  #                     
  #                     Cunningham <- f_Cunningham(rad_species)
  #                     Diffusivity <- f_Diffusivity(Matrix = "air", Temp, DynViscAirStandard, rad_species, Cunningham)
  #                     
  #                     SchmidtNumber <- DynViscAirStandard / (rhoMatrix * Diffusivity) # rhoMatrix to be converted to RhoWater or RhoAir
  #                     # alpha.surf = depends on vegetation type, see e.g. LOTEUR ref guide table A.1
  #                     if (to.Matrix == "water") {
  #                       Brown <- SchmidtNumber^-(1/2)
  #                     } else 
  #                     {
  #                       Brown <- SchmidtNumber^(-2/3) 
  #                       }
  #                   
  #                     rad_RainDrop <- f_RadRain(RAINrate, FRACtwet)
  #                     
  #                     Cunningham.cw <- f_Cunningham(rad_RainDrop)
  #                     rho_water <- 998
  #                     Settvel.Particle.cw <-  ((2*rad_RainDrop)^2 * (rho_water - from.rhoMatrix)* #rhoMatrix to be converted to RhoWater and RhoAir
  #                                                  constants::syms$gn*Cunningham.cw)/(18*DynViscAirStandard) 
  #                     
  #                     Relax.Particle.a <- ((rho_species-from.rhoMatrix)*(2*rad_species)^2*f_Cunningham(rad_species))/(18*DynViscAirStandard)
  #                     
  #                     StN <- (2*Relax.Particle.a*(Settvel.Particle.cw-SettlingVelocity))/(2*rad_RainDrop)
  #                     
  #                     AREAFRACveg <- 0.01
  #                     LargeVegRadius <- 5e-4
  #                     VegHair <- 1e-5
  #                     
  #                     Intercept <- 0.3*(AREAFRACveg * (rad_species/(rad_species+VegHair))+(1-AREAFRACveg)*(rad_species/(rad_species+LargeVegRadius)))
  #                     
  #                     # StokesNumberVeg <- fStokesNumber_rain(rad_particle, rho_species, from.rho, from.visc, RAINrate, FRACtwet, Cunningham.cw, g) # to be update in future!
  #                     R1 <- exp(-StN^0.5) # R1 = correction factor representing the fraction of particles that stick to the surface
  #                     epsilon = 3 # epsilon is empirical constant set to 3
  #                     beta.a = 2 # constant set to 2
  #                     
  #                     # in SB4N alpha.surf is 0.8!!!
  #                     if (to.Matrix == "water") {
  #                       Impaction <- 10^(-3/StN)
  #                     } else if (to.Matrix == "soil") {
  #                       Impaction <- (StN / (StN + 0.8))^beta.a
  #                     }
  #                     
  #                     if (to.Matrix == "water") {
  #                       SurfResist <- 1 /(FricVel * (Brown + Impaction))
  #                     } else {
  #                       SurfResist <- 1 / (FricVel * (Brown + # Collection efficiency for Brownian diffusion
  #                                                       Intercept + # Collection efficiency for interception
  #                                                       Impaction))
  #                     }
  #                     sea_AERO <- 135
  #                     if (to.SubCompartName == "sea" ){
  #                       DRYDEPvelocity <- 1 /(sea_AERO + SurfResist) + SettlingVelocity
  #                     } else
  #                     {
  #                     
  #                     # Collection efficiency for impaction
  #                     DRYDEPvelocity <- 1 / (AEROresist + SurfResist) + SettlingVelocity
  #                     }
  #                     
  #                     # Currently implemented in SimpleBox for P species:
  #                     if (SpeciesName == "Attached") {
  #                       DRYDEPvelocity <- AEROSOLdeprate # AEROSOLdeprate constant given in xls version of SB4
  #                     }
  #                     return((DRYDEPvelocity * to.Area) / from.Volume)
  #                     
  #                   },
  #                   NA)
  #            
  #          } else {
  #            return(NA)
  #          }
  #        },
  #        "FALSE" = {
  #          # Existing function logic when test is FALSE
  #          if (SpeciesName %in% c("Nanoparticle", "Aggregated", "Attached")) {
  #            if (anyNA(c(AEROresist, DynViscAirStandard, rhoMatrix, rho_species, to.alpha.surf))) {
  #              return(NA)
  #            }
  #            switch(SubCompartName,
  #                   "air" = {
  #                     
  #                     Cunningham <- f_Cunningham(rad_species)
  #                     Diffusivity <- f_Diffusivity(Matrix = "air", Temp, DynViscAirStandard, rad_species, Cunningham)
  #                     
  #                     SchmidtNumber <- DynViscAirStandard / (rhoMatrix * Diffusivity) # rhoMatrix to be converted to RhoWater or RhoAir
  #                     # gamma.surf = depends on vegetation type, see e.g. LOTEUR ref guide table A.1
  #                     Brown <- SchmidtNumber^(-to.gamma.surf)
  #                     
  #                     StN <- ifelse(to.ColRad == 0 | is.na(to.ColRad), # StokesNumber following ref guide LOTEUR v2.0 2016
  #                                   (SettlingVelocity * FricVel) / DynViscAirStandard / rhoMatrix, # for smooth surfaces (water)
  #                                   (SettlingVelocity * FricVel) / (constants::syms$gn * to.ColRad)      # for vegetated surfaces (soil)
  #                     )
  #                     Intercept <- ifelse(to.ColRad == 0 | is.na(to.ColRad),
  #                                         0, # for smooth surfaces
  #                                         0.5 * (rad_species / to.ColRad)^2 # LOTEUR (eq 5.14) in ref guide.
  #                     )
  #                     
  #                     # StokesNumberVeg <- fStokesNumber_rain(rad_particle, rho_species, from.rho, from.visc, RAINrate, FRACtwet, Cunningham.cw, g) # to be update in future!
  #                     R1 <- exp(-StN^0.5) # R1 = correction factor representing the fraction of particles that stick to the surface
  #                     epsilon = 3 # epsilon is empirical constant set to 3
  #                     beta.a = 2 # constant set to 2
  #                     
  #                     # in SB4N alpha.surf is 0.8!!!
  #                     Impaction <- (StN / (StN + to.alpha.surf))^beta.a # LOTUS EUROS ref guide (eq. 5.13)  
  #                     
  #                     SurfResist <- 1 / (epsilon * FricVel * (Brown + # Collection efficiency for Brownian diffusion
  #                                                               Intercept + # Collection efficiency for interception
  #                                                               Impaction) * R1) # Collection efficiency for impaction
  #                     DRYDEPvelocity <- 1 / (AEROresist + SurfResist) + SettlingVelocity
  #                     
  #                     # Currently implemented in SimpleBox for P species:
  #                     if (SpeciesName == "Attached") {
  #                       DRYDEPvelocity <- AEROSOLdeprate # AEROSOLdeprate constant given in xls version of SB4
  #                     }
  #                     return((DRYDEPvelocity * to.Area) / from.Volume)
  #                     
  #                   },
  #                   NA)
  #            
  #          } else {
  #            return(NA)
  #          }
  #        }
  # )
}
