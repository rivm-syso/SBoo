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
rho_species <- function (SpeciesName, SubCompartName,ScaleName,
                         RhoS, RadS, 
                         RadCOL, RadCP, RhoCOL, RhoCP,
                         RhoNuc, RadNuc, 
                         NumConcNuc, NumConcAcc,
                         Shortest_side, Intermediate_side,
                         Longest_side, Shape){
  if (((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & (SubCompartName == "freshwatersediment" | 
                                                             SubCompartName == "lakesediment" |
                                                             SubCompartName == "lake" |
                                                             SubCompartName == "river" |
                                                             SubCompartName == "agriculturalsoil"|
                                                             SubCompartName == "othersoil")) | 
      (ScaleName %in% c("Regional", "Continental")) & (SubCompartName == "deepocean" )) {
    return(NA)
  }
  if (is.na(Shape) || is.null(Shape)) {
    Shape <- "Default"
  }
  
  if (Shape == "Sphere" | Shape == "Default") {
    if(is.na(RadS)||is.null(RadS)){RadS = Shortest_side/2}
    switch(SpeciesName,
           "Nanoparticle" = return (RhoS),
           "Aggregated" = {
             if(SubCompartName == "air") {
               SingleMass <- ((NumConcNuc*(RhoNuc*fVol(RadNuc)+RhoS*fVol(RadS)))+(NumConcAcc*(RhoCOL*(fVol(RadCOL))+RhoS*fVol(RadS)))) /
                 (NumConcNuc+NumConcAcc)
               SingleVol <- ((NumConcNuc*(fVol(RadS)+((fVol(RadNuc)))))+(NumConcAcc*(fVol(RadS)+(fVol(RadCOL))))) /
                 (NumConcNuc+NumConcAcc)
               return(SingleMass/SingleVol)
             } else {
               SingleMass <- RhoS*fVol(RadS) + RhoCOL*fVol(RadCOL) # No more DF with 1/3
               SingleVol <-  fVol(RadS) + fVol(RadCOL)
               return(SingleMass/SingleVol)
             }
           },
           "Attached" = {
             SingleMass <- RhoS*fVol(RadS) + RhoCP*fVol(RadCP) # No more DF with 1/3
             SingleVol <- fVol(RadCP) + fVol(RadS)
             return(SingleMass/SingleVol)
           },
           return(NA)
    )
    
    
  } else if (Shape == "Ellipsoid" | Shape == "Cube" | Shape == "Box" | 
             Shape == "Film" | Shape == "Cylindric - circular" | 
             Shape == "Cylindric - elliptic") {
    switch(SpeciesName,
           "Nanoparticle" = return(RhoS),
           "Aggregated" = {
             if(SubCompartName == "air") {
               SingleMass <- 
                 ((NumConcNuc * 
                     (RhoS * fVol(Shape = Shape,
                                  Longest_side = Longest_side,
                                  Intermediate_side = Intermediate_side,
                                  Shortest_side = Shortest_side) +
                        RhoNuc * fVol(RadNuc))) +
                    (NumConcAcc * 
                       (RhoS * fVol(Shape = Shape,
                                    Longest_side = Longest_side,
                                    Intermediate_side = Intermediate_side,
                                    Shortest_side = Shortest_side) +
                          RhoCOL * fVol(RadCOL)))) / 
                 (NumConcNuc + NumConcAcc)
               
               SingleVol <- ((NumConcNuc * (fVol(
                 Shape = Shape,
                 Longest_side = Longest_side,
                 Intermediate_side = Intermediate_side,
                 Shortest_side = Shortest_side
               ) +
                 fVol(RadNuc))) +
                 (NumConcAcc * (fVol(
                   Shape = Shape,
                   Longest_side = Longest_side,
                   Intermediate_side = Intermediate_side,
                   Shortest_side = Shortest_side
                 ) +
                   fVol(RadCOL)))) / (NumConcNuc + NumConcAcc)
               return(SingleMass/SingleVol)
             } else {
               SingleMass <- RhoS*fVol(Shape = Shape,
                                       Longest_side = Longest_side,
                                       Intermediate_side = Intermediate_side,
                                       Shortest_side = Shortest_side)  + RhoCOL*fVol(RadCOL)
               SingleVol <- fVol(Shape = Shape,
                                 Longest_side = Longest_side,
                                 Intermediate_side = Intermediate_side,
                                 Shortest_side = Shortest_side) + fVol(RadCOL)
               return(SingleMass/SingleVol)
             }
           },
           "Attached" = {
             SingleMass <- RhoS*fVol(Shape = Shape,
                                     Longest_side = Longest_side,
                                     Intermediate_side = Intermediate_side,
                                     Shortest_side = Shortest_side) + RhoCP*fVol(RadCP)
             SingleVol <- fVol(RadCP) + fVol(Shape = Shape,
                                             Longest_side = Longest_side,
                                             Intermediate_side = Intermediate_side,
                                             Shortest_side = Shortest_side)
             return(SingleMass/SingleVol)
           },
           return(NA)
    )
  } else {
    return("Invalid Shape! Please choose from Sphere, Ellipsoid, Cube, Box, Cylindric - circular, or Cylindric - elliptic.")
  }
  
  
}
