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
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return rho_species, approach to calculate density of species in air or water [kg.m-3]
#' @export
rho_species <- function (SpeciesName, SubCompartName,ScaleName,
                         RhoS, RadS, 
                         RadCOL, RadCP, RhoCOL, RhoCP,
                         RhoNuc, RadNuc, 
                         NumConcNuc, NumConcAcc,
                         Shortest_side, Intermediate_side,
                         Longest_side, Shape, Regional_and_Continental_deepocean,
                         parent){
  if (!Shape %in% c("Ellipsoid", "Cube", "Box", "Film", "Cylindric - circular", "Cylindric - elliptic", "Fiber", "Sphere", "Default")) {
    stop("Invalid Shape! Please choose from Sphere, Ellipsoid, Cube, Box, Cylindric - circular, Fiber, or Cylindric - elliptic.")
  }
  
  out <- ScaleName |>
    tidyr::expand_grid(SubCompartName, SpeciesName) |>
    dplyr::left_join(RadCOL, by=c("SubCompart")) |>
    dplyr::left_join(RhoCOL, by=c("SubCompart")) |>
    dplyr::left_join(RadCP, by=c("SubCompart")) |>
    dplyr::left_join(RhoCP, by=c("SubCompart")) |>
    dplyr::left_join(NumConcNuc, by=c("Scale")) |>
    dplyr::left_join(NumConcAcc, by=c("Scale")) |>
    parent$states$clipStates() |>
    dplyr::mutate(
      RhoS = RhoS, RadS = RadS,RhoNuc = RhoNuc,RadNuc =RadNuc, Shape = Shape, Shortest_side=Shortest_side,Intermediate_side=Intermediate_side, Longest_side=Longest_side,
      # Some conditions to set defaults for input values
      Shape = dplyr::case_when(
        is.na(Shape) | is.null(Shape) ~ 'Default',
        TRUE ~ Shape
      ),
      RadS = dplyr::case_when(
        Shape %in% c("Sphere", "Default") & (is.na(RadS) | is.null(RadS)) ~ Shortest_side/2,
        TRUE ~ RadS
      )
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      Vol_RadS = fVol(RadS),
      Vol_RadNuc = fVol(RadNuc),
      Vol_RadCOL = fVol(RadCOL),
      Vol_RadCP = fVol(RadCP),
      Vol_shape = fVol(rad_particle=RadS, Shape = Shape,Longest_side = Longest_side,Intermediate_side = Intermediate_side,Shortest_side = Shortest_side)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      rho_species = dplyr::case_when(
        SpeciesName == 'Nanoparticle' ~ RhoS, 
        (ScaleName %in% c("Tropic", "Moderate", "Arctic") & SubCompartName %in% c("freshwatersediment", "lakesediment", "lake", "river", "agriculturalsoil", "othersoil")) |
          (ScaleName %in% c("Regional", "Continental") & SubCompartName == "deepocean" & (isFALSE(Regional_and_Continental_deepocean) | is.na(Regional_and_Continental_deepocean) | Regional_and_Continental_deepocean == "FALSE")) ~
          NA_real_,
        Shape %in% c("Sphere", "Default") & SpeciesName == 'Aggregated' & SubCompart == 'air' ~
          (((NumConcNuc*(RhoNuc*Vol_RadNuc+RhoS*Vol_RadS))+(NumConcAcc*(RhoCOL*(Vol_RadCOL)+RhoS*Vol_RadS))) /
          (NumConcNuc+NumConcAcc)) / (((NumConcNuc*(Vol_RadS+((Vol_RadNuc))))+(NumConcAcc*(Vol_RadS+(Vol_RadCOL)))) /
          (NumConcNuc+NumConcAcc)),
        
        Shape %in% c("Sphere", "Default") & SpeciesName == 'Aggregated' & SubCompart != 'air' ~
          (RhoS*Vol_RadS + RhoCOL*Vol_RadCOL) / (Vol_RadS + Vol_RadCOL),
        
        Shape %in% c("Sphere", "Default") & SpeciesName == 'Attached' ~ 
          (RhoS*Vol_RadS + RhoCP*Vol_RadCP) / (Vol_RadCP + Vol_RadS),
        
        Shape %in% c("Ellipsoid", "Cube", "Box", "Film", "Cylindric - circular", "Cylindric - elliptic", "Fiber") & SpeciesName == 'Aggregated' & SubCompart == 'air'~
          (((NumConcNuc * (RhoS * Vol_shape + RhoNuc * Vol_RadNuc)) + (NumConcAcc * (RhoS * Vol_shape + RhoCOL * Vol_RadCOL))) / (NumConcNuc + NumConcAcc)) /
          (((NumConcNuc * (Vol_shape +Vol_RadNuc)) +(NumConcAcc * (Vol_shape +Vol_RadCOL))) / (NumConcNuc + NumConcAcc)),
        
        Shape %in% c("Ellipsoid", "Cube", "Box", "Film", "Cylindric - circular", "Cylindric - elliptic", "Fiber") & SpeciesName == 'Aggregated' & SubCompart != 'air'~
          (RhoS*Vol_shape  + RhoCOL*Vol_RadCOL) / (Vol_shape + Vol_RadCOL),
        
        Shape %in% c("Ellipsoid", "Cube", "Box", "Film", "Cylindric - circular", "Cylindric - elliptic", "Fiber") & SpeciesName == 'Attached'~
          (RhoS*Vol_shape + RhoCP*Vol_RadCP) / (Vol_RadCP + Vol_shape),
        
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(rho_species)) |>
    dplyr::select(Scale, SubCompart, Species, rho_species) |>
    dplyr::arrange(Scale, SubCompart, Species)
  
  return(data.frame(out))
  
  # if (((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & (SubCompartName == "freshwatersediment" | 
  #                                                            SubCompartName == "lakesediment" |
  #                                                            SubCompartName == "lake" |
  #                                                            SubCompartName == "river" |
  #                                                            SubCompartName == "agriculturalsoil"|
  #                                                            SubCompartName == "othersoil")) | 
  #     ((ScaleName %in% c("Regional", "Continental")) & (SubCompartName == "deepocean" ) &&
  #      (isFALSE(Regional_and_Continental_deepocean) || is.na(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "FALSE"))) {
  #   return(NA)
  # }
  # if (is.na(Shape) || is.null(Shape)) {
  #   Shape <- "Default"
  # }
  # 
  # if (Shape == "Sphere" | Shape == "Default") {
  #   if(is.na(RadS)||is.null(RadS)){RadS = Shortest_side/2}
  #   switch(SpeciesName,
  #          "Nanoparticle" = return (RhoS),
  #          "Aggregated" = {
  #            if(SubCompartName == "air") {
  #              SingleMass <- ((NumConcNuc*(RhoNuc*fVol(RadNuc)+RhoS*fVol(RadS)))+(NumConcAcc*(RhoCOL*(fVol(RadCOL))+RhoS*fVol(RadS)))) /
  #                (NumConcNuc+NumConcAcc)
  #              SingleVol <- ((NumConcNuc*(fVol(RadS)+((fVol(RadNuc)))))+(NumConcAcc*(fVol(RadS)+(fVol(RadCOL))))) /
  #                (NumConcNuc+NumConcAcc)
  #              return(SingleMass/SingleVol)
  #            } else {
  #              SingleMass <- RhoS*fVol(RadS) + RhoCOL*fVol(RadCOL) # No more DF with 1/3
  #              SingleVol <-  fVol(RadS) + fVol(RadCOL)
  #              return(SingleMass/SingleVol)
  #            }
  #          },
  #          "Attached" = {
  #            SingleMass <- RhoS*fVol(RadS) + RhoCP*fVol(RadCP) # No more DF with 1/3
  #            SingleVol <- fVol(RadCP) + fVol(RadS)
  #            return(SingleMass/SingleVol)
  #          },
  #          return(NA)
  #   )
  #   
  #   
  # } else if (Shape == "Ellipsoid" | Shape == "Cube" | Shape == "Box" | 
  #            Shape == "Film" | Shape == "Cylindric - circular" | 
  #            Shape == "Cylindric - elliptic" | Shape == "Fiber") {
  #   switch(SpeciesName,
  #          "Nanoparticle" = return(RhoS),
  #          "Aggregated" = {
  #            if(SubCompartName == "air") {
  #              SingleMass <- 
  #                ((NumConcNuc * 
  #                    (RhoS * fVol(rad_particle=RadS,
  #                                 Shape = Shape,
  #                                 Longest_side = Longest_side,
  #                                 Intermediate_side = Intermediate_side,
  #                                 Shortest_side = Shortest_side) +
  #                       RhoNuc * fVol(RadNuc))) +
  #                   (NumConcAcc * 
  #                      (RhoS * fVol(rad_particle=RadS,
  #                                   Shape = Shape,
  #                                   Longest_side = Longest_side,
  #                                   Intermediate_side = Intermediate_side,
  #                                   Shortest_side = Shortest_side) +
  #                         RhoCOL * fVol(RadCOL)))) / 
  #                (NumConcNuc + NumConcAcc)
  #              
  #              SingleVol <- ((NumConcNuc * (fVol(
  #                rad_particle=RadS,
  #                Shape = Shape,
  #                Longest_side = Longest_side,
  #                Intermediate_side = Intermediate_side,
  #                Shortest_side = Shortest_side
  #              ) +
  #                fVol(RadNuc))) +
  #                (NumConcAcc * (fVol(
  #                  rad_particle=RadS,
  #                  Shape = Shape,
  #                  Longest_side = Longest_side,
  #                  Intermediate_side = Intermediate_side,
  #                  Shortest_side = Shortest_side
  #                ) +
  #                  fVol(RadCOL)))) / (NumConcNuc + NumConcAcc)
  #              return(SingleMass/SingleVol)
  #            } else {
  #              SingleMass <- RhoS*fVol(rad_particle=RadS,
  #                                      Shape = Shape,
  #                                      Longest_side = Longest_side,
  #                                      Intermediate_side = Intermediate_side,
  #                                      Shortest_side = Shortest_side)  + RhoCOL*fVol(RadCOL)
  #              SingleVol <- fVol(rad_particle=RadS,
  #                                Shape = Shape,
  #                                Longest_side = Longest_side,
  #                                Intermediate_side = Intermediate_side,
  #                                Shortest_side = Shortest_side) + fVol(RadCOL)
  #              return(SingleMass/SingleVol)
  #            }
  #          },
  #          "Attached" = {
  #            SingleMass <- RhoS*fVol(rad_particle=RadS,
  #                                    Shape = Shape,
  #                                    Longest_side = Longest_side,
  #                                    Intermediate_side = Intermediate_side,
  #                                    Shortest_side = Shortest_side) + RhoCP*fVol(RadCP)
  #            SingleVol <- fVol(RadCP) + fVol(rad_particle=RadS,
  #                                            Shape = Shape,
  #                                            Longest_side = Longest_side,
  #                                            Intermediate_side = Intermediate_side,
  #                                            Shortest_side = Shortest_side)
  #            return(SingleMass/SingleVol)
  #          },
  #          return(NA)
  #   )
  # } else {
  #   return("Invalid Shape! Please choose from Sphere, Ellipsoid, Cube, Box, Cylindric - circular, Fiber, or Cylindric - elliptic.")
  # }
  
  
}
