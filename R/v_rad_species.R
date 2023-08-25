#' @title Radius of species
#' @name rad_species
#' @description Calculate radius of heteroagglomerates [m]
#' @param SpeciesName description
#' @param SubCompartName description
#' @param NaturalRad natural particle radius all but small in air [m]
#' @param RadS nanoparticle radius [m]
#' @param Df fractal dimension of combined heteroagglomerate [-]
#' @param RadNuc Nucleation mode aerosol particle radius [m]
#' @param RadAcc Accumulation mode aerosol particle radius [m]
#' @param NumConcNuc Number concentration of Nucleation mode aerosol particles [#/m3]
#' @param NumConcAcc Number concentration of Accumulation mode aerosol particles [#/m3]
#' @param ... dimension(s, Species, Scales, SubCompart) of the
#' @return rad_species, Approach to calculate the radius of small heteroagglomerates in air/water [m]
#' @export
rad_species <- function(SpeciesName,SubCompartName,NaturalRad,RadNuc,RadAcc,RadS,NumConcNuc,NumConcAcc, Df) {
  if(anyNA(c(RadNuc,RadAcc,NumConcNuc,NumConcAcc))) return(NA)

  if (SpeciesName == "Aggregated" & SubCompartName == "air"){ # 2 types of aerosols (Acc & Nuc) combining to become a heteroagglomerate with and ENP.
    SingleVol <- ((NumConcNuc*(fVol(RadS)+fVol(RadNuc)))+(NumConcAcc*(fVol(RadS)+fVol(RadAcc))))/(NumConcNuc+NumConcAcc)
    rad_particle <- (SingleVol/((4/3)*pi))^(Df)
    return(rad_particle)
  } else {
    if (SpeciesName == "Nanoparticle") {
      return(RadS)
    } else
    {
      return((NaturalRad^3 + RadS^3)^(Df))
    }
  }
}
