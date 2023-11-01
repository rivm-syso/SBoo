#' @title Heteroagglomeration with aerosols
#' @name k_HeteroAgglomeration.a
#' @description Calculation of the first order rate constant (s-1) for heteroagglomeration of ENPs with large (>0.45 um) natural particles
#' @param to.alpha Attachment Efficiency of ENPs with other particulates [-]
#' @param NumConc_Otherparticle Mass concentration of other particulates [kg.m-3]
#' @param RadS Radius of nanoparticle [m]
#' @param RhoS Density of nanoparticle [kg.m-3]
#' @param radius_Otherparticle Radius of natural particle [m]
#' 
#' 
#' @param rho_Otherparticle Density (specific weight) of natural particle [kg/m3]
#' @param rhoMatrix Density of fluid matrix (water) [kg/m3]

#' @param Temp Temperature of compartment [K]
#' @param viscosity Dyanic viscosity of fluid matrix []
#' @return k_HeteroAgglomeration, the rate constant for 1rst order process: heteroagglomeration [s-1]
#' @export
#'
k_HeteroAgglomeration.a <- function(rad_species, RadOther, rho_species, RhoOther, Temp,
                                    Diffusivity, 
                                    NumConcNuc, NumConcAcc, NumConcCP){

  ThermVel <- function(Temp, Radius, Rho){
   ((8*KBOLTZ*Temp)/(pi*fVol(Radius)*Rho))^0.5
  }
  
  f.Fuchs <- function(Diff.P1,
                      Diff.P2,
                      Rad.P1,
                      Rad.P2,
                      Thermvel.P1,
                      Thermvel.P2){
    (1+(4*(Diff.P1+Diff.P2))/((Rad.P1+Rad.P2)*(Thermvel.P1+Thermvel.P2^2)^0.5))^-1
  }
  
  ThermVelSpecies <- ThermVel(Temp, rad_species, rho_species)
  radLarge <-  FetchOtherDim("rad_species", Species = "Large")
  if (rad_species >= radLarge) { #else: for Small particle, it's more complicated
    DiffOther <- fDiffusivity(Temp = Temp, visc = DynViscAirStandard, rad_particle = RadOther)
    ThermVelOther <- ThermVel(Temp, RadOther, RhoOther)
    Fuchs <- f.Fuchs (Diffusivity, DiffOther,
                      rad_species, RadOther,
                      ThermVelSpecies, ThermVelOther)
    return(Fuchs*(4*PI()*(RadOther+rad_species)*(Diff.Other+Diffusivity))*NumConcCP)
  } else { #rad_species < radLarge ??Solid also ?!
    # for Acc
    DiffAcc <- fDiffusivity(Temp = Temp, visc = DynViscAirStandard, rad_particle = RadAcc)
    ThermVelAcc <- ThermVel(Temp, RadAcc, RhoAcc)
    FuchsAcc <- f.Fuchs (Diffusivity, DiffAcc,
                      rad_species, RadAcc,
                      ThermVelSpecies, ThermVelAcc)
    k_Acc <- FuchsAcc*(4*PI()*(RadAcc+rad_species)*(DiffNuc+Diffusivity))*NumConcAcc
    # for Nuc
    DiffNuc <- fDiffusivity(Temp = Temp, visc = DynViscAirStandard, rad_particle = RadNuc)
    ThermVelNuc <- ThermVel(Temp, RadNuc, RhoNuc)
    FuchsNuc <- f.Fuchs (Diffusivity, DiffNuc,
                         rad_species, RadNuc,
                         ThermVelSpecies, ThermVelNuc)
    k_Nuc <- Fuchs*(4*PI()*(RadNuc+rad_species)*(DiffNuc+Diffusivity))*NumConcNuc
    return(k_Acc + k_Nuc)
  }
}
