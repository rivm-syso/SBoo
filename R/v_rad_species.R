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
#' @return rad_species, Approach to calculate the radius of small heteroagglomerates in air/water [m]
#' @export
rad_species <- function(SpeciesName, SubCompartName,ScaleName,
                        RadCOL, RadCP, RadNuc,
                        RadS, NumConcNuc, NumConcAcc,
                        Shortest_side, Intermediate_side,
                        Longest_side, Shape) {
  
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
    switch(tolower(SpeciesName),
      "nanoparticle" = return(RadS),
      "aggregated" = {
        if (tolower(SubCompartName) == "air") {
          SingleVol <- ((NumConcNuc * (fVol(RadS) + fVol(RadNuc))) + (NumConcAcc * (fVol(RadS) + fVol(RadCOL)))) / (NumConcNuc + NumConcAcc)
          rad_particle <- (SingleVol / ((4 / 3) * pi))^(1 / 3)
          return(rad_particle)
        } else {
          SingleVol <- fVol(RadS) + fVol(RadCOL)
          rad_particle <- (SingleVol / ((4 / 3) * pi))^(1 / 3)
          return(rad_particle)
        }
      },
      "attached" = {
        SingleVol <- fVol(RadS) + fVol(RadCP)
        rad_particle <- (SingleVol / ((4 / 3) * pi))^(1 / 3)
        return(rad_particle)
      },
      return(NA)
    )
  } else if (Shape == "Ellipsoid") {
    switch(tolower(SpeciesName),
      "nanoparticle" = return(Shortest_side),
      "aggregated" = {
        if (tolower(SubCompartName) == "air") {
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
          Shortest_side <- SingleVol / ((1 / 6) * pi * Intermediate_side * Longest_side)
          return(Shortest_side / 2)
        } else {
          SingleVol <- fVol(
            Shape = Shape,
            Longest_side = Longest_side,
            Intermediate_side = Intermediate_side,
            Shortest_side = Shortest_side
          ) + fVol(RadCOL)
          Shortest_side <- SingleVol / ((1 / 6) * pi * Intermediate_side * Longest_side)
          return(rad_particle)
        }
      },
      "attached" = {
        SingleVol <- fVol(
          Shape = Shape,
          Longest_side = Longest_side,
          Intermediate_side = Intermediate_side,
          Shortest_side = Shortest_side
        ) + fVol(RadCP)
        Shortest_side <- SingleVol / ((1 / 6) * pi * Intermediate_side * Longest_side)
        return(rad_particle)
      },
      return(NA)
    )
  } else if (Shape == "Cube" | Shape == "Box" | Shape == "Film") {
    switch(tolower(SpeciesName),
      "nanoparticle" = return(Shortest_side),
      "aggregated" = {
        if (tolower(SubCompartName) == "air") {
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
          Shortest_side <- SingleVol / (Intermediate_side * Longest_side)
          return(Shortest_side / 2)
        } else {
          SingleVol <- fVol(
            Shape = Shape,
            Longest_side = Longest_side,
            Intermediate_side = Intermediate_side,
            Shortest_side = Shortest_side
          ) + fVol(RadCOL)
          Shortest_side <- SingleVol / (Intermediate_side * Longest_side)
          return(Shortest_side / 2)
        }
      },
      "attached" = {
        SingleVol <- fVol(
          Shape = Shape,
          Longest_side = Longest_side,
          Intermediate_side = Intermediate_side,
          Shortest_side = Shortest_side
        ) + fVol(RadCP)
        Shortest_side <- SingleVol / (Intermediate_side * Longest_side)
        return(Shortest_side / 2)
      },
      return(NA)
    )
  } else if (Shape == "Cylindric - circular") {
    switch(tolower(SpeciesName),
      "nanoparticle" = return(Shortest_side),
      "aggregated" = {
        if (tolower(SubCompartName) == "air") {
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
          Shortest_side <- (SingleVol / (Longest_side * pi))^(1 / 2)
          return(Shortest_side / 2)
        } else {
          SingleVol <- fVol(
            Shape = Shape,
            Longest_side = Longest_side,
            Intermediate_side = Intermediate_side,
            Shortest_side = Shortest_side
          ) + fVol(RadCOL)
          Shortest_side <- (SingleVol / (Longest_side * pi))^(1 / 2)
          return(Shortest_side / 2)
        }
      },
      "attached" = {
        SingleVol <- fVol(
          Shape = Shape,
          Longest_side = Longest_side,
          Intermediate_side = Intermediate_side,
          Shortest_side = Shortest_side
        ) + fVol(RadCP)
        Shortest_side <- (SingleVol / (Longest_side * pi))^(1 / 2)
        return(Shortest_side / 2)
      },
      return(NA)
    )
  } else if (Shape == "Cylindric - elliptic") {
    radius_major <- Intermediate_side / 2
    radius_minor <- Shortest_side / 2
    height <- Longest_side
    volume <- pi * radius_major * radius_minor * height
    return(volume)

    switch(tolower(SpeciesName),
      "nanoparticle" = return(Shortest_side),
      "aggregated" = {
        if (tolower(SubCompartName) == "air") {
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
          rad_minor_particle <- (SingleVol / (Longest_side / 2 * Intermediate_side / 2 * pi))
          return(rad_minor_particle)
        } else {
          SingleVol <- fVol(
            Shape = Shape,
            Longest_side = Longest_side,
            Intermediate_side = Intermediate_side,
            Shortest_side = Shortest_side
          ) + fVol(RadCOL)
          rad_minor_particle <- (SingleVol / (Longest_side / 2 * Intermediate_side / 2 * pi))
          return(rad_minor_particle)
        }
      },
      "attached" = {
        SingleVol <- fVol(
          Shape = Shape,
          Longest_side = Longest_side,
          Intermediate_side = Intermediate_side,
          Shortest_side = Shortest_side
        ) + fVol(RadCP)
        rad_minor_particle <- (SingleVol / (Longest_side / 2 * Intermediate_side / 2 * pi))
        return(rad_minor_particle)
      },
      return(NA)
    )
  } else {
    return(stop("Invalid Shape! Please choose from Sphere, Ellipsoid, Cube, Box, Cylindric - circular, or Cylindric - elliptic."))
  }
}
