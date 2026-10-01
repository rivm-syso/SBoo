#' @title Radius of species
#' @name rad_species
#' @description Calculate radius of heteroagglomerates [m]
#' @param SpeciesName species considered
#' @param SubCompartName subcompartment considered
#' @param NaturalRad natural particle radius all but small in air [m]
#' @param RadS nanoparticle radius [m]
#' @param Shortest_side Shortest side of the particle (nanomaterial or microplastics) this side will be adjusted based on the natural particle [m]
#' @param Intermediate_side Intermediate side of the particle (nanomaterial or microplastics) [m]
#' @param Longest_side Longest side of the particle (nanomaterial or microplastics) [m]
#' @param Shape Shape of the particle
#' @param RadNuc Nucleation mode aerosol particle radius [m]
#' @param RadCOL Accumulation mode aerosol particle radius [m]
#' @param RadCP coarse particulate mode aerosol particle radius [m]
#' @param NumConcNuc Number concentration of Nucleation mode aerosol particles [#/m3]
#' @param NumConcAcc Number concentration of Accumulation mode aerosol particles [#/m3]
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return rad_species, Approach to calculate the radius of small heteroagglomerates in air/water [m]
#' @export
rad_species <- function(SpeciesName, SubCompartName,ScaleName,
                        RadCOL, RadCP, RadNuc,
                        RadS, NumConcNuc, NumConcAcc,
                        Shortest_side, Intermediate_side,
                        Longest_side, Shape, Regional_and_Continental_deepocean,
                        parent) {
  
  # First check and do warnings/stops
  if (Shape %in% c('Ellipsoid', 'Cube', 'Box', 'Film', 'Cylindric - elliptic') & (is.na(Intermediate_side) ||
                              is.null(Intermediate_side) ||
                              is.na(Longest_side) ||
                              is.null(Longest_side) ||
                              is.null(Shortest_side) ||
                              is.null(Shortest_side))) {
    stop(paste("For", Shape, "shape a side dimension is missing (e.g. Intermediate_side)"))
  }
  if (Shape %in% c("Cylindric - circular", "Fiber") & (is.na(Longest_side) ||
                                                          is.null(Longest_side) ||
                                                          is.null(Shortest_side) ||
                                                          is.null(Shortest_side))) {
    stop(paste("For", Shape, "shape a side dimension is missing (e.g. Longest_side"))
  }

  out <- ScaleName |>
    tidyr::expand_grid(SubCompartName, SpeciesName) |>
    parent$states$clipStates() |>
    dplyr::left_join(RadCOL, by=c("SubCompart")) |>
    dplyr::left_join(RadCP, by=c("SubCompart")) |>
    dplyr::left_join(NumConcNuc, by=c("Scale")) |>
    dplyr::left_join(NumConcAcc, by=c("Scale")) |>
    dplyr::mutate(
      # Some input constants
      Shortest_side = Shortest_side,
      Intermediate_side = Intermediate_side,
      Longest_side = Longest_side,
      Shape = Shape,
      RadNuc = RadNuc,
      RadS = RadS,

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
    # Calculate volumes rowwise
    dplyr::rowwise() |>
    dplyr::mutate(
      Vol_RadS = fVol(RadS),
      Vol_RadNuc = fVol(RadNuc),
      Vol_RadCOL = fVol(RadCOL),
      Vol_RadCP = fVol(RadCP),
      Vol_shape = fVol(Shape = Shape,Longest_side = Longest_side,Intermediate_side = Intermediate_side,Shortest_side = Shortest_side)
    ) |>
    dplyr::ungroup() |>
    # Determine rad species
    dplyr::mutate(
      rad_species = dplyr::case_when(
        # 1st : if any compartment in global or deepocean in regional and continental == NA (any shape)
        (ScaleName %in% c("Tropic", "Moderate", "Arctic") & SubCompartName %in% c("freshwatersediment", "lakesediment", "lake", "river", "agriculturalsoil", "othersoil")) |
          (ScaleName %in% c("Regional", "Continental") & SubCompartName == "deepocean" & (isFALSE(Regional_and_Continental_deepocean) | is.na(Regional_and_Continental_deepocean) | Regional_and_Continental_deepocean == "FALSE")) ~
          NA,

        # 2nd : if its a Sphere or Default and Nanoparticle return RadS
        Shape %in% c("Sphere", "Default") & SpeciesName == 'Nanoparticle' ~ RadS,

        # 3rd : if its a Sphere or Default and aggregate and air return calculated rad particle
        Shape %in% c("Sphere", "Default") & SpeciesName == 'Aggregated' & SubCompart == 'air' ~
          (((NumConcNuc * (Vol_RadS + Vol_RadNuc)) + (NumConcAcc * (Vol_RadS + Vol_RadCOL))) / (NumConcNuc + NumConcAcc) / ((4 / 3) * pi))^(1 / 3),

        # 4th : if its a Sphere or Default and aggregate and NOT air return calculated rad particle
        Shape %in% c("Sphere", "Default") & SpeciesName == 'Aggregated' & SubCompart != 'air' ~
          ((Vol_RadS + Vol_RadCOL) / ((4 / 3) * pi))^(1 / 3),

        # 5th : if its a Sphere or Default and attached and ANY subcompart return calculated rad particle
        Shape %in% c("Sphere", "Default") & SpeciesName == 'Attached' ~
          ((Vol_RadS + Vol_RadCP) / ((4 / 3) * pi))^(1 / 3),
        
        # 6th : if sphere but anything else
        Shape %in% c("Sphere", "Default") & !SpeciesName  %in% c("Attached", "Aggregated", "Nanoparticle") ~
          NA,

        # 7th : if its an ellipsoid and nanoparticle return shortest side
        Shape == 'Ellipsoid'  & SpeciesName == 'Nanoparticle' ~ Shortest_side,

        # 8th : if its an ellipsoid and aggregate and subcompart is air return calculated rad_s
        Shape == 'Ellipsoid'  & SpeciesName == 'Aggregated' & SubCompart == 'air' ~
          ((((NumConcNuc * (Vol_shape +Vol_RadNuc)) + (NumConcAcc * (Vol_shape +Vol_RadCOL))) / (NumConcNuc + NumConcAcc)) / ((1 / 6) * pi * Intermediate_side * Longest_side)) / 2,

        # 9th : if its an ellipsoid and aggregate and subcompart is NOT air return calculated rad_s
        Shape == 'Ellipsoid'  & SpeciesName == 'Aggregated' & SubCompart != 'air' ~
          ((Vol_shape + Vol_RadCOL)/ ((1 / 6) * pi * Intermediate_side * Longest_side)) / 2,
        
        # 10th : if its an ellipsoid and attached and subcompart return calc rad_s
        Shape == 'Ellipsoid'  & SpeciesName == 'Attached' ~
          ((Vol_shape + Vol_RadCOL)/ ((1 / 6) * pi * Intermediate_side * Longest_side)) / 2,
        
        # 11th : if ellipsoid and anything else
        Shape == 'Ellipsoid' & !SpeciesName  %in% c("Attached", "Aggregated", "Nanoparticle") ~
          NA,
  
        # 12th : if cube, box or film and nanoparticle
        Shape %in% c("Cube", "Box", "Film") & SpeciesName == 'Nanoparticle' ~ Shortest_side,
        
        # 13th : if cube, box or film and aggregate and subcompart is air
        Shape %in% c("Cube", "Box", "Film") & SpeciesName == 'Aggregated' & SubCompart == 'air' ~
          ((((NumConcNuc * (Vol_shape + Vol_RadNuc)) + (NumConcAcc * (Vol_shape + Vol_RadCOL))) / (NumConcNuc + NumConcAcc)) / (Intermediate_side * Longest_side)) /2,

        # 14th : if cube, box or film and aggregate and subcompart is NOT air return calculated rad_s
        Shape %in% c("Cube", "Box", "Film") & SpeciesName == 'Aggregated' & SubCompart != 'air' ~
        ((Vol_shape + Vol_RadCOL)/ (Intermediate_side * Longest_side)) / 2,
  
        # 15th : if cube, box or film and attached and subcompart is any return calculated rad_s
        Shape %in% c("Cube", "Box", "Film") & SpeciesName == 'Attached' ~ 
          ((Vol_shape + Vol_RadCP) / (Intermediate_side * Longest_side)) / 2,
        
        # 16th : if Cube box or film all else return  NA
        Shape %in% c("Cube", "Box", "Film") & !SpeciesName  %in% c("Attached", "Aggregated", "Nanoparticle") ~
          NA,
  
        # 16th : if cylindric or fiber and nanoparticle return shortest side
        Shape %in% c("Cylindric - circular", "Fiber") & SpeciesName == 'Nanoparticle' ~ Shortest_side,
  
        # 13th : if "Cylindric - circular", "Fiber" and aggregate and subcompart is air
        Shape %in% c("Cylindric - circular", "Fiber") & SpeciesName == 'Aggregated' & SubCompart == 'air'~ 
          ((((NumConcNuc * (Vol_shape + Vol_RadNuc)) + (NumConcAcc * (Vol_shape + Vol_RadCOL))) / (NumConcNuc + NumConcAcc)) / (Longest_side * pi))^(1 / 2),
  
        # 14th : if "Cylindric - circular", "Fiber" and aggregate and subcompart is NOT air
        Shape %in% c("Cylindric - circular", "Fiber") & SpeciesName == 'Aggregated' & SubCompart != 'air'~ 
          ((Vol_shape + Vol_RadCOL)/ (Longest_side * pi))^(1/2),
        
        # 15th : if "Cylindric - circular", "Fiber" and Attached
        Shape %in% c("Cylindric - circular", "Fiber") & SpeciesName == 'Attached'~ 
          ((Vol_shape + Vol_RadCOL)/ (Longest_side * pi))^(1/2),
  
        # 16th : if Cylindric - elliptic and nanoparticle return shortest side
        Shape == "Cylindric - elliptic" & SpeciesName == 'Nanoparticle' ~ Shortest_side,
      
        # 17th : if cylindric - elliptic and aggregate and subcompart is air return calc rad s
        Shape == "Cylindric - elliptic" & SpeciesName == 'Aggregated' & SubCompart == 'air' ~
          (((((NumConcNuc * (Vol_shape + Vol_RadNuc)) + (NumConcAcc * (Vol_shape + Vol_RadCOL))) / (NumConcNuc + NumConcAcc))) / (Longest_side / 2 * Intermediate_side / 2 * pi)),
  
        # 18th : if cylindric - elliptic and aggregate and subcompart is NOT air return calc rad s
        Shape == "Cylindric - elliptic" & SpeciesName == 'Aggregated' & SubCompart != 'air' ~
          ((Vol_shape + Vol_RadCOL)/ (Longest_side / 2 * Intermediate_side / 2 * pi)),
  
        # 19th : if cylindric - elliptic and Attached
        Shape %in% c("Cylindric - circular", "Fiber") & SpeciesName == 'Attached' ~ 
          ((Vol_shape + Vol_RadCP)/ (Longest_side / 2 * Intermediate_side / 2 * pi)),
        
        TRUE ~ NA
    )
    ) |>
    dplyr::filter(!is.na(rad_species)) |>
    dplyr::select(Scale, SubCompart, Species, rad_species) |>
    dplyr::arrange(Scale, SubCompart, Species)
  
  return(data.frame(out))
# 
#   if (
#     (ScaleName %in% c("Tropic", "Moderate", "Arctic") & 
#      SubCompartName %in% c("freshwatersediment", "lakesediment", "lake", "river", "agriculturalsoil", "othersoil")
#     ) | 
#     (
#       ScaleName %in% c("Regional", "Continental") & 
#       SubCompartName == "deepocean" & 
#       (isFALSE(Regional_and_Continental_deepocean) | is.na(Regional_and_Continental_deepocean) | Regional_and_Continental_deepocean == "FALSE")
#     )
#   ) {
#     return(NA)
#   }
#   
#   if (is.na(Shape) || is.null(Shape)) {
#     Shape <- "Default"
#   }
# 
# 
#   if (Shape == "Sphere" | Shape == "Default") {
#     if(is.na(RadS)||is.null(RadS)){RadS = Shortest_side/2}
#     switch(tolower(SpeciesName),
#       "nanoparticle" = return(RadS),
#       "aggregated" = {
#         if (tolower(SubCompartName) == "air") {
#           SingleVol <- ((NumConcNuc * (fVol(RadS) + fVol(RadNuc))) + (NumConcAcc * (fVol(RadS) + fVol(RadCOL)))) / (NumConcNuc + NumConcAcc)
#           rad_particle <- (SingleVol / ((4 / 3) * pi))^(1 / 3)
#           return(rad_particle)
#         } else {
#           SingleVol <- fVol(RadS) + fVol(RadCOL)
#           rad_particle <- (SingleVol / ((4 / 3) * pi))^(1 / 3)
#           return(rad_particle)
#         }
#       },
#       "attached" = {
#         SingleVol <- fVol(RadS) + fVol(RadCP)
#         rad_particle <- (SingleVol / ((4 / 3) * pi))^(1 / 3)
#         return(rad_particle)
#       },
#       return(NA)
#     )
#   } else if (Shape == "Ellipsoid") {
#     
#     if(is.na(Intermediate_side) || 
#        is.null(Intermediate_side) ||
#        is.na(Longest_side) || 
#        is.null(Longest_side) ||
#        is.null(Shortest_side) ||
#        is.null(Shortest_side)) stop(paste("For", Shape, "shape a side dimension is missing (e.g. Intermediate_side)"))
#     
#     switch(tolower(SpeciesName),
#       "nanoparticle" = return(Shortest_side),
#       "aggregated" = {
#         if (tolower(SubCompartName) == "air") {
#           SingleVol <- ((NumConcNuc * (fVol(
#             # rad_particle=RadS, For ellipsoid shortest side should be given, not rads
#             Shape = Shape,
#             Longest_side = Longest_side,
#             Intermediate_side = Intermediate_side,
#             Shortest_side = Shortest_side
#           ) +
#             fVol(RadNuc))) +
#             (NumConcAcc * (fVol(
#               # rad_particle=RadS,
#               Shape = Shape,
#               Longest_side = Longest_side,
#               Intermediate_side = Intermediate_side,
#               Shortest_side = Shortest_side
#             ) +
#               fVol(RadCOL)))) / (NumConcNuc + NumConcAcc)
#           Shortest_side <- SingleVol / ((1 / 6) * pi * Intermediate_side * Longest_side)
#           return(Shortest_side / 2)
#         } else {
#           SingleVol <- fVol(
#             # rad_particle=RadS,
#             Shape = Shape,
#             Longest_side = Longest_side,
#             Intermediate_side = Intermediate_side,
#             Shortest_side = Shortest_side
#           ) + fVol(RadCOL)
#           Shortest_side <- SingleVol / ((1 / 6) * pi * Intermediate_side * Longest_side)
#           return(Shortest_side)
#         }
#       },
#       "attached" = {
#         SingleVol <- fVol(
#           # rad_particle=RadS,
#           Shape = Shape,
#           Longest_side = Longest_side,
#           Intermediate_side = Intermediate_side,
#           Shortest_side = Shortest_side
#         ) + fVol(RadCP)
#         Shortest_side <- SingleVol / ((1 / 6) * pi * Intermediate_side * Longest_side)
#         return(Shortest_side)
#       },
#       return(NA)
#     )
#   } else if (Shape == "Cube" | Shape == "Box" | Shape == "Film") {
#     
#     if(is.na(Intermediate_side) || 
#        is.null(Intermediate_side) ||
#        is.na(Longest_side) || 
#        is.null(Longest_side) ||
#        is.null(Shortest_side) ||
#        is.null(Shortest_side)) stop(paste("For",Shape,"shape a side dimension is missing (e.g. Intermediate_side)"))
#     
#     switch(tolower(SpeciesName),
#       "nanoparticle" = return(Shortest_side),
#       "aggregated" = {
#         if (tolower(SubCompartName) == "air") {
#           SingleVol <- ((NumConcNuc * (fVol(
#             # rad_particle=RadS,
#             Shape = Shape,
#             Longest_side = Longest_side,
#             Intermediate_side = Intermediate_side,
#             Shortest_side = Shortest_side
#           ) +
#             fVol(RadNuc))) +
#             (NumConcAcc * (fVol(
#               # rad_particle=RadS,
#               Shape = Shape,
#               Longest_side = Longest_side,
#               Intermediate_side = Intermediate_side,
#               Shortest_side = Shortest_side
#             ) +
#               fVol(RadCOL)))) / (NumConcNuc + NumConcAcc)
#           Shortest_side <- SingleVol / (Intermediate_side * Longest_side)
#           return(Shortest_side / 2)
#         } else {
#           SingleVol <- fVol(
#             # rad_particle=RadS,
#             Shape = Shape,
#             Longest_side = Longest_side,
#             Intermediate_side = Intermediate_side,
#             Shortest_side = Shortest_side
#           ) + fVol(RadCOL)
#           Shortest_side <- SingleVol / (Intermediate_side * Longest_side)
#           return(Shortest_side / 2)
#         }
#       },
#       "attached" = {
#         SingleVol <- fVol(
#           # rad_particle=RadS,
#           Shape = Shape,
#           Longest_side = Longest_side,
#           Intermediate_side = Intermediate_side,
#           Shortest_side = Shortest_side
#         ) + fVol(RadCP)
#         Shortest_side <- SingleVol / (Intermediate_side * Longest_side)
#         return(Shortest_side / 2)
#       },
#       return(NA)
#     )
#   } else if (Shape == "Cylindric - circular" | Shape == "Fiber") {
#     if(is.na(Longest_side) || 
#        is.null(Longest_side) ||
#        is.null(Shortest_side) ||
#        is.null(Shortest_side)) stop(paste("For",Shape,"shape a side dimension is missing (e.g. Longest_side)"))
#     
#     switch(tolower(SpeciesName),
#       "nanoparticle" = return(Shortest_side),
#       "aggregated" = {
#         if (tolower(SubCompartName) == "air") {
#           SingleVol <- ((NumConcNuc * (fVol(
#             # rad_particle=RadS,
#             Shape = Shape,
#             Longest_side = Longest_side,
#             Shortest_side = Shortest_side
#           ) +
#             fVol(RadNuc))) +
#             (NumConcAcc * (fVol(
#               # rad_particle=RadS,
#               Shape = Shape,
#               Longest_side = Longest_side,
#               Shortest_side = Shortest_side
#             ) +
#               fVol(RadCOL)))) / (NumConcNuc + NumConcAcc)
#           Shortest_side <- (SingleVol / (Longest_side * pi))^(1 / 2)
#           return(Shortest_side / 2)
#         } else {
#           SingleVol <- fVol(
#             # rad_particle=RadS,
#             Shape = Shape,
#             Longest_side = Longest_side,
#             Shortest_side = Shortest_side
#           ) + fVol(RadCOL)
#           Shortest_side <- (SingleVol / (Longest_side * pi))^(1 / 2)
#           return(Shortest_side / 2)
#         }
#       },
#       "attached" = {
#         SingleVol <- fVol(
#           # rad_particle=RadS,
#           Shape = Shape,
#           Longest_side = Longest_side,
#           Shortest_side = Shortest_side
#         ) + fVol(RadCP)
#         Shortest_side <- (SingleVol / (Longest_side * pi))^(1 / 2)
#         return(Shortest_side / 2)
#       },
#       return(NA)
#     )
#   } else if (Shape == "Cylindric - elliptic") {
#     if(is.na(Intermediate_side) || 
#        is.null(Intermediate_side) ||
#        is.na(Longest_side) || 
#        is.null(Longest_side) ||
#        is.null(Shortest_side) ||
#        is.null(Shortest_side)) stop(paste("For",Shape,"shape a side dimension is missing (e.g. Intermediate_side)"))
# 
#     switch(tolower(SpeciesName),
#       "nanoparticle" = return(Shortest_side),
#       "aggregated" = {
#         if (tolower(SubCompartName) == "air") {
#           SingleVol <- ((NumConcNuc * (fVol(
#             # rad_particle=RadS,
#             Shape = Shape,
#             Longest_side = Longest_side,
#             Intermediate_side = Intermediate_side,
#             Shortest_side = Shortest_side
#           ) +
#             fVol(RadNuc))) +
#             (NumConcAcc * (fVol(
#               # rad_particle=RadS,
#               Shape = Shape,
#               Longest_side = Longest_side,
#               Intermediate_side = Intermediate_side,
#               Shortest_side = Shortest_side
#             ) +
#               fVol(RadCOL)))) / (NumConcNuc + NumConcAcc)
#           rad_minor_particle <- (SingleVol / (Longest_side / 2 * Intermediate_side / 2 * pi))
#           
#           if(is_na(rad_minor_particle)) stop("Result from v_rad_species is NA")
#           
#           return(rad_minor_particle)
#         } else {
#           SingleVol <- fVol(
#             # rad_particle=RadS,
#             Shape = Shape,
#             Longest_side = Longest_side,
#             Intermediate_side = Intermediate_side,
#             Shortest_side = Shortest_side
#           ) + fVol(RadCOL)
#           rad_minor_particle <- (SingleVol / (Longest_side / 2 * Intermediate_side / 2 * pi))
#           return(rad_minor_particle)
#         }
#       },
#       "attached" = {
#         SingleVol <- fVol(
#           # rad_particle=RadS,
#           Shape = Shape,
#           Longest_side = Longest_side,
#           Intermediate_side = Intermediate_side,
#           Shortest_side = Shortest_side
#         ) + fVol(RadCP)
#         rad_minor_particle <- (SingleVol / (Longest_side / 2 * Intermediate_side / 2 * pi))
#         return(rad_minor_particle)
#       },
#       return(NA)
#     )
#   } else {
#     return(stop("Invalid Shape! Please choose from Sphere, Ellipsoid, Cube, Box, Cylindric - circular, Fiber, or Cylindric - elliptic."))
#   }
}
