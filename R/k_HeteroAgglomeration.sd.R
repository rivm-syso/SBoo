#' @title Formation of heteroagglomerate (P) of ENPs (S) and solid grains in soil or sediment
#' @name k.HeteroAgglomeration.sd
#' @description Calculation of the first order rate constant (s-1) for heteroagglomeration of ENPs with large (>0.45 um) natural particles
#' @param alpha Attachment Efficiency of ENPs with other particulates [-]
#' @param rad_species Radius (of nanoparticle) [m]
#' @param rho_species Density (of nanoparticle) [kg.m-3]
#' @param hamakerSP.w Hamaker constant heteroagglomerate (ENP, water, SiO2) [J]
#' @param RadCP Radius of natural particle associated with the "to" species [m]
#' @param Udarcy Darcu velocity in soil or sediment [m.s-1]
#' @param FRACs Fraction solid particles in soil or sediment [-]
#' @param Temp Temperature of Scale [K]
#' @return k.HeteroAgglomeration, the rate constant for 1rst order process: heteroagglomeration [s-1]
#' @export
k_HeteroAgglomeration.sd <- function(Udarcy, rad_species, rho_species, 
                                     RadCP, 
                                     to.FRACs, to.alpha,
                                     Temp, hamakerSP.w, DynViscWaterStandard,
                                     Matrix){

  DiffS.w <- f_Diffusivity(Matrix=Matrix, 
                           Temp, DynVisc=DynViscWaterStandard, 
                           rad_species=rad_species)

  rhoWater <- 998
  Por <- 1-to.FRACs
  GammPDF <- (1-Por)^(1/3)
  
  ASPDF <- (2*(1-GammPDF^5))/(2-3*GammPDF+3*GammPDF^5-2*GammPDF^6)
  aspectratioSFP <- rad_species/RadCP
  PecletNumberFP <- (Udarcy*2*RadCP)/(DiffS.w)
  vdWaalsNumberSFP <- hamakerSP.w/(KBOLTZ*Temp)
  
  BrownSFP <- 2.4*ASPDF^(1/3)*aspectratioSFP^(-0.081)*PecletNumberFP^-0.715*vdWaalsNumberSFP^0.053
  
  InterceptSFP <- 0.55*aspectratioSFP^1.55*PecletNumberFP^-0.125*vdWaalsNumberSFP^0.125
  
  GravNumberS <- (2*rad_species^2*(rho_species-rhoWater)*g)/(9*DynViscWaterStandard*Udarcy)
  GravSFP <- 0.22*aspectratioSFP^-0.24*GravNumberS^1.11*vdWaalsNumberSFP^0.053
  
  fTotalSFP <- BrownSFP+InterceptSFP+GravSFP

  Filter <- (3/2)*(1-Por)/(2*RadCP*Por)

  K_het.sd <- Filter*Udarcy*fTotalSFP*to.alpha 
  
  return(K_het.sd)
}
