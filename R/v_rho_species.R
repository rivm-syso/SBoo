#' @title Rho (density) of species
#' @name rho_species
#' @param SpeciesName Dimension
#' @param SubCompartName Dimension
#' @param RadS nanoparticle radius [m]
#' @param RhoS nanoparticle density [kg/m3]
#' @param Df fractal dimension of combined heteroagglomerate [-]
#' @param RadNuc Nucleation mode aerosol particle radius [m]
#' @param RadAcc Accumulation mode aerosol particle radius [m]
#' @param RhoNuc Nucleation mode aerosol particle  [kg/m3]
#' @param RhoAcc Accumulation mode aerosol particle  [kg/m3]
#' @param NumConcNuc Number concentration of Nucleation mode aerosol particles [#/m3]
#' @param NumConcAcc Number concentration of Accumulation mode aerosol particles [#/m3]
#' @return rho_species, approach to calculate density of species in air or water [kg.m-3]
#' @export
rho_species <- function (SpeciesName, SubCompartName,
                         RhoS, RadS, 
                         NaturalRad, NaturalRho,
                         RhoNuc, RadNuc,
                         RhoAcc, RadAcc,
                         NumConcNuc, NumConcAcc, Df){
  
  if (SpeciesName == "Nanoparticle"){ 
    return (RhoS)
  }
  if (SpeciesName == "Molecular"){
    return(NA)
  }
  if (SpeciesName == "Aggregated"  & SubCompartName == "air"){ 
    #count weighted average for nucleation mode and Aitken accumulation mode 
    SingleMass <- ((NumConcNuc*(RhoNuc*fVol(RadNuc)+RhoS*fVol(RadS)))+(NumConcAcc*(RhoAcc*(fVol(RadAcc))+RhoS*fVol(RadS)))) /
      (NumConcNuc+NumConcAcc)
    SingleVol <- ((NumConcNuc*(fVol(RadS)+((fVol(RadNuc)))))+(NumConcAcc*(fVol(RadS)+(fVol(RadAcc))))) /
      (NumConcNuc+NumConcAcc)
  }
  else {
    SingleMass <- RhoS*fVol(RadS) + NaturalRho*fVol(NaturalRad) #Requires update for DF is not 1/3, include matrix mass
    SingleVol <- fVol((NaturalRad^3 + RadS^3)^(Df)) 
  }
  SingleMass/SingleVol
}
