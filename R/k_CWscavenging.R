#' @title Cloud water scavenging coefficients as function of particle size and density
#' @name k_CWscavenging
#' @description Process describing the transport from dry air to rain, based on Wang et al. (2010) https://doi.org/10.5194/acp-10-5685-2010
#' @param rho_species particle density [kg m-3]
#' @param rad_species particle radius [m]
#' @param rad_Raindrop radius of the raindrop [m]
#' @param to.rhoMatrix density of cloud water [kg m-3]
#' @param rhoMatrix density of air [kg m-3]
#' @param DynViscAirStandard dynamic viscosity of air [kg m-3 s-1]
#' @param DynViscWaterStandard dynamic viscosity of (cloud) water [kg m-3 s-1]
#' @param Temp temperature [K]
#' @param RAINrate Yearly rain rate [m/s]
#' @param FRACtwet Time fraction of wet episodes [-]
#' @param tdry duration of dry episodes [s]
#' @param twet duration of wet episodes [s]
#' @param COLLECTeff Aerosol collection efficiency for raindrops (Wang et al., 2010) [-]
#' @param SettlingVelocity settling velocity for particle species in air [m/s]
#' @param Matrix compartment type 
#' @param SpeciesName name of the considered species 
#' @param SubCompartName considered subcompartment 
#' @param VertDistance depth of compartment [m]
#' @param CritStokesNumb.cw Critical stokes number, based on Wang et al. (2010)
#' @param ReyNumber.cw Reynolds number 
#' @param SchmidtNumber.Particle.a Denotes the ratio of momentum diffusivity and mass diffusivity, based on Wang et al. (2010)
#' @return k_CWscavenging
#' @export
k_CWscavenging <- function(RAINrate, FRACtwet, tdry, twet, COLLECTeff,
                           rad_species, rho_species, rhoMatrix, 
                           DynViscAirStandard, DynViscWaterStandard, SettlingVelocity,
                           Temp, Matrix, SpeciesName, VertDistance,
                           SubCompartName, ScaleName, parent){
  
  
  # Custom join func
  join_suffix <- function(x, y, by=NULL) {
    dplyr::left_join(x, y, by = by, suffix = c("_from", "_to"))
  }

  # Names of dataframe
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_CWscavenging")))
  
  # Build grid
  args <- c("ScaleName", "SubCompartName", "Matrix", "SpeciesName", "RAINrate", "FRACtwet", "tdry", "twet", "COLLECTeff", "SettlingVelocity", "rad_species", "rho_species", "rhoMatrix", "Temp", "VertDistance")
  dfs <- mget(args, envir = environment())
  grid <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  
  out <- parent$FromDataAndTo("k_CWscavenging") |>
    # From
    dplyr::left_join(grid, by=c("Scale" = "Scale", "from.SubCompart" = "SubCompart", "Species" = "Species")) |>
    
    # To met functie om suffix in te stellen
    join_suffix(Matrix, by=c("to.SubCompart" = "SubCompart")) |>
    join_suffix(rhoMatrix, by=c("Matrix_to" = "Matrix")) |>

    dplyr::mutate(
      rad_RainDrop = f_RadRain(RAINrate, FRACtwet),
      Cunningham.cw = f_Cunningham(rad_RainDrop),
      Cunningham_rad_species = f_Cunningham(rad_species),
      Settvel.Particle.cw =  2*(rad_RainDrop^2 * (rhoMatrix_to - rhoMatrix_from)* #rhoMatrix to be converted to RhoWater and RhoAir
                                   constants::syms$gn*Cunningham.cw)/(9*DynViscAirStandard),
      Relax.Particle.a = ((rho_species-rhoMatrix_from)*(2*rad_species)^2*Cunningham_rad_species)/(18*DynViscAirStandard),
      StokesNumber.Particle.a = (2*Relax.Particle.a*(Settvel.Particle.cw-SettlingVelocity))/(2*rad_RainDrop),
      ReyNumber.cw = ((2*rad_RainDrop)*Settvel.Particle.cw*rhoMatrix_from)/(2*DynViscAirStandard),
      SchmidtNumber.Particle.a = DynViscAirStandard/(rhoMatrix_from*f_Diffusivity(Matrix_from,Temp,DynViscAirStandard,rad_species,Cunningham = Cunningham_rad_species)),
      CritStokesNumb.cw = ((1.2+(1/12)*log(1+ReyNumber.cw))/(1+log(1+ReyNumber.cw))),
      Grav.a.cw = dplyr::case_when(
        StokesNumber.Particle.a > CritStokesNumb.cw ~ ((StokesNumber.Particle.a-CritStokesNumb.cw)/(StokesNumber.Particle.a-CritStokesNumb.cw+2/3))^(3/2),
        TRUE ~ 0
      ),
      Intercept.a.cw = 4*(rad_species/rad_RainDrop)*((DynViscAirStandard/DynViscWaterStandard)+(1+2*ReyNumber.cw^0.5*(rad_species/rad_RainDrop))),
      Brown.a.cw = (4/(ReyNumber.cw*SchmidtNumber.Particle.a))*(1+0.4*(ReyNumber.cw^0.5*SchmidtNumber.Particle.a^(1/3))+0.16*(ReyNumber.cw^0.5*SchmidtNumber.Particle.a^0.5)),
      Total = Brown.a.cw + Intercept.a.cw + Grav.a.cw,
      
      k_CWscavenging = dplyr::case_when(
        is.na(StokesNumber.Particle.a) ~ NA,
        from.SubCompart != "air" ~ NA,
        SpeciesName %in% c("Nanoparticle","Aggregated") ~ (3/2)*(Total*RAINrate)/(2*rad_RainDrop),
        SpeciesName == "Attached" ~ ((tdry+twet)/twet*RAINrate*COLLECTeff)/VertDistance,
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(k_CWscavenging)) |>
    dplyr::select(names, k_CWscavenging) |>
    dplyr::arrange(dplyr::across(names))

  return(data.frame(out))
  
  # if(SubCompartName != "air") return(NA)
  # 
  # if(SpeciesName %in% c("Nanoparticle","Aggregated")){
  # # variables for calculation of 3 types of collection mechanisms (Gravitational, Intercept, Brownian)
  # 
  # rad_RainDrop <- f_RadRain(RAINrate, FRACtwet)
  # 
  # Cunningham.cw <- f_Cunningham(rad_RainDrop)
  # Settvel.Particle.cw <-  2*(rad_RainDrop^2 * (to.rhoMatrix - from.rhoMatrix)* #rhoMatrix to be converted to RhoWater and RhoAir
  #                              constants::syms$gn*Cunningham.cw)/(9*DynViscAirStandard) 
  # 
  # Relax.Particle.a <- ((rho_species-from.rhoMatrix)*(2*rad_species)^2*f_Cunningham(rad_species))/(18*DynViscAirStandard)
  # 
  # StokesNumber.Particle.a <- (2*Relax.Particle.a*(Settvel.Particle.cw-SettlingVelocity))/(2*rad_RainDrop)
  # 
  # ReyNumber.cw <- ((2*rad_RainDrop)*Settvel.Particle.cw*from.rhoMatrix)/(2*DynViscAirStandard)
  # 
  # SchmidtNumber.Particle.a <- DynViscAirStandard/(from.rhoMatrix*f_Diffusivity(Matrix, 
  #                                                                       Temp, 
  #                                                                       DynViscAirStandard, 
  #                                                                       rad_species,  
  #                                                                       Cunningham = f_Cunningham(rad_species)))
  # 
  # CritStokesNumb.cw <- ((1.2+(1/12)*log(1+ReyNumber.cw))/(1+log(1+ReyNumber.cw)))
  # 
  # if (is.na(StokesNumber.Particle.a)) return(NA)
  # 
  # # calculation of graviational collection efficiency
  # if(StokesNumber.Particle.a > CritStokesNumb.cw){
  #   Grav.a.cw <- ((StokesNumber.Particle.a-CritStokesNumb.cw)/(StokesNumber.Particle.a-CritStokesNumb.cw+2/3))^(3/2)
  # } else Grav.a.cw <- 0
  # 
  # # calculation of Interception collection efficiency
  # Intercept.a.cw <- 4*(rad_species/rad_RainDrop)*
  #   ((DynViscAirStandard/DynViscWaterStandard)+(1+2*ReyNumber.cw^0.5*(rad_species/rad_RainDrop)))
  # 
  # # calculation of Brownian collection efficiency
  # 
  # Brown.a.cw <- (4/(ReyNumber.cw*SchmidtNumber.Particle.a))*
  #   (1+0.4*(ReyNumber.cw^0.5*SchmidtNumber.Particle.a^(1/3))+0.16*(ReyNumber.cw^0.5*SchmidtNumber.Particle.a^0.5))
  # #from f_Brown??  ((2*getConst("r")*Temp)/(3*viscosity))*(from.radius+radOther )^2/(from.radius*radOther)
  # 
  # Total <- Brown.a.cw + Intercept.a.cw + Grav.a.cw
  # 
  # # RAINrate.wet <- RAINrate/FRACtwet
  # 
  # (3/2)*(Total*RAINrate)/(2*rad_RainDrop)
  # } else if (SpeciesName == "Attached") {
  #   ((tdry+twet)/twet*RAINrate*COLLECTeff)/VertDistance
  # } else return(NA)
}


