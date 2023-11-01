#' @title Cloud water scavenging coefficients as function of particle size and density
#' @name k_CWscavenging
#' @description a description follows .. Joris? based on Wang, X.; Zhang, L.;
#'   Moran, M.D. Uncertainty assessment of current size-resolved
#'   parameterizations for below-cloud particle scavenging by rain. Atmos. Chem.
#'   Phys. Discuss. 2010, 10, 2503–2548
#' @param rho_species particle density []
#' @param rad_species particle radius []
#' @param to.rhoMatrix density of cloud water []
#' @param rhoMatrix density of air []
#' @param DynViscAirStandard dynamic viscosity of air []
#' @param DynViscWaterStandard dynamic viscosity of (cloud) water []
#' @param Temp tempreature [oC]
#' @param 
#' @return k_CWscavenging
#' @export
k_CWscavenging <- function(RAINrate, FRACtwet, tdry, twet, COLLECTeff,
                           rad_species, rho_species, to.rhoMatrix, from.rhoMatrix, 
                           DynViscAirStandard, DynViscWaterStandard, from.SettlingVelocity,
                           Temp, Matrix, SpeciesName, VertDistance,
                           SubCompartName){
  
  if(SubCompartName != "air") return(NA)
  
  if(SpeciesName %in% c("Nanoparticle","Aggregated")){
  # variables for calculation of 3 types of collection mechanisms (Gravitational, Intercept, Brownian)
  
  rad_RainDrop <- f_RadRain(RAINrate, FRACtwet)
  
  Cunningham.cw <- f_Cunningham(rad_RainDrop)
  Settvel.Particle.cw <-  2*(rad_RainDrop^2 * (to.rhoMatrix - from.rhoMatrix)*
                               constants::syms$gn*Cunningham.cw)/(9*DynViscAirStandard) 
  
  Relax.Particle.a <- ((rho_species-from.rhoMatrix)*(2*rad_species)^2*f_Cunningham(rad_species))/(18*DynViscAirStandard)
  
  StokesNumber.Particle.a  = (2*Relax.Particle.a*(Settvel.Particle.cw-from.SettlingVelocity))/(2*rad_RainDrop)
  
  ReyNumber.cw =((2*rad_RainDrop)*Settvel.Particle.cw*from.rhoMatrix)/(2*DynViscAirStandard)
  
  SchmidtNumber.Particle.a = DynViscAirStandard/(from.rhoMatrix*f_Diffusivity(Matrix, 
                                                                        Temp, 
                                                                        DynViscAirStandard, 
                                                                        rad_species,  
                                                                        Cunningham = f_Cunningham(rad_species)))
  
  CritStokesNumb.cw = ((1.2+(1/12)*log(1+ReyNumber.cw))/(1+log(1+ReyNumber.cw)))
  
  # calculation of graviational collection efficiency
  if(StokesNumber.Particle.a>CritStokesNumb.cw){
    Grav.a.cw = ((StokesNumber.Particle.a-CritStokesNumb.cw)/(StokesNumber.Particle.a-CritStokesNumb.cw+2/3))^(3/2)
  } else Grav.a.cw = 0
  
  # calculation of Interception collection efficiency
  Intercept.a.cw <- 4*(rad_species/rad_RainDrop)*
    ((DynViscAirStandard/DynViscWaterStandard)+(1+2*ReyNumber.cw^0.5*(rad_species/rad_RainDrop)))
  
  # calculation of Brownian collection efficiency
  
  Brown.a.cw <- (4/(ReyNumber.cw*SchmidtNumber.Particle.a))*
    (1+0.4*(ReyNumber.cw^0.5*SchmidtNumber.Particle.a^(1/3))+0.16*(ReyNumber.cw^0.5*SchmidtNumber.Particle.a^0.5))
  #from f_Brown??  ((2*getConst("r")*Temp)/(3*viscosity))*(from.radius+radOther )^2/(from.radius*radOther)
  
  Total <- Brown.a.cw + Intercept.a.cw + Grav.a.cw
  
  # RAINrate.wet <- RAINrate/FRACtwet
  
  (3/2)*(Total*RAINrate)/(2*rad_RainDrop)
  } else if (SpeciesName == "Attached") {
    ((tdry+twet)/twet*RAINrate*COLLECTeff)/VertDistance
  } else return(NA)
}


