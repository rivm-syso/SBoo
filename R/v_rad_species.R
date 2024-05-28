#' @title Radius of species
#' @name rad_species
#' @description Calculate radius of heteroagglomerates [m]
#' @param SpeciesName description
#' @param SubCompartName description
#' @param NaturalRad natural particle radius all but small in air [m]
#' @param RadS nanoparticle radius [m]
#' @param Df fractal dimension of combined heteroagglomerate [-]
#' @param RadNuc Nucleation mode aerosol particle radius [m]
#' @param RadCOL Accumulation mode aerosol particle radiusy [m]
#' @param NumConcNuc Number concentration of Nucleation mode aerosol particles [#/m3]
#' @param NumConcAcc Number concentration of Accumulation mode aerosol particles [#/m3]
#' @param ... dimension(s, Species, Scales, SubCompart) of the
#' @return rad_species, Approach to calculate the radius of small heteroagglomerates in air/water [m]
#' @export
rad_species <- function(SpeciesName,SubCompartName,
                        RadCOL,RadCP,RadNuc,
                        RadS,NumConcNuc,NumConcAcc, Df) {
  # if(anyNA(c(RadNuc,RadCOL,NumConcNuc,NumConcAcc))) return(NA)
  # browser()
  switch (tolower(SpeciesName),
          "nanoparticle" = return(RadS),
          "aggregated" = {
            if (tolower(SubCompartName) == "air") {
              SingleVol <- ((NumConcNuc*(fVol(RadS)+fVol(RadNuc)))+(NumConcAcc*(fVol(RadS)+fVol(RadCOL))))/(NumConcNuc+NumConcAcc)
              rad_particle <- (SingleVol/((4/3)*pi))^(Df)
              return(rad_particle)
              #return((rad_particle^3 + RadS^3)^(Df))
            } else return((RadCOL^3 + RadS^3)^(Df))
          },
          "attached" = return((RadCP^3 + RadS^3)^(Df)),
          return(NA)
  )
}
