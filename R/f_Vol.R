#' @title Calculate volume of sphere or other shape
#' @name fVol
#' @description Calculate the volume of spherical particle using the radius in m3, shapes based on Gov4Nano 
#' @param rad_particle Radius of particle [m]
#' @param Shape Particle shape defined by user, used for nano and plastic [-]
#' @param Longest_side Longest side length as defined by user [m]
#' @param Intermediate_side Intermediate side length as defined by user [m]
#' @param Shortest_side Shortest side length as defined by user
#' @return fVol [m3]
#' @export

fVol <- function(rad_particle, #option to use input as radius
                 Shape = NULL, 
                 Longest_side = NULL, # same characteristic as radius above, but given as diameter or longest side of other shape
                 Intermediate_side = NULL, 
                 Shortest_side = NULL){
  if (is.na(Shape) || is.null(Shape)){
    Shape <- "Default"
  }
  
  # Check if Shortest side is NA or NULL and assign default values if so
  if ( is.na(Shortest_side) || is.null(Shortest_side) ) {
    Shortest_side <- rad_particle * 2
  }
  # Check if any of Intermediate or Longest sides is NA or NULL and assign default values if so
  if (is.na(Intermediate_side) || is.null(Intermediate_side) ||is.na(Longest_side) || is.null(Longest_side)) {
    Intermediate_side <- rad_particle * 2 #maybe 0.75 or build in shape functions
    Longest_side <- rad_particle * 2
  }
  
  
  if (Shape == "Sphere" | Shape == "Default") {
    radius <- Shortest_side / 2
    volume <- (4/3) * pi * radius^3
    return(volume)
  } else if (Shape == "Ellipsoid") {
    volume <- (1/6) * pi * Longest_side * Intermediate_side * Shortest_side # with sides being full length, not half.
    return(volume)
  } else if (Shape == "Cube" | Shape == "Box" | Shape == "Film") {
    #Longest_side <- sqrt(3)*Longest_side
    #Intermediate_side <-sqrt(2)*Longest_side
    volume <- Longest_side * Intermediate_side * Shortest_side
    return(volume)
  } else if (Shape == "Cylindric - circular" | Shape == "Fiber") {
    
    radius <- Shortest_side / 2
    height <- Longest_side
    volume <- pi * radius^2 * height
    return(volume)
  } else if (Shape == "Cylindric - elliptic") {
    radius_major <- Intermediate_side  / 2
    radius_minor <- Shortest_side  / 2
    height <- Longest_side
    volume <- pi * radius_major * radius_minor * height
    return(volume)
  } else {
    return("Invalid Shape! Please choose from Sphere, Ellipsoid, Cube, Box, Film, Fiber, Cylindric - circular, or Cylindric - elliptic.")
  }
}

