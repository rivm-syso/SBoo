#'@title SettVellNat
#'@name SettVellNat
#'@description Settling velocity for natural suspended particles
#'@param NaturalRad
#'@param NaturalRho
#'@param rhoMatrix
#'@param DynVisc
#'@export
SettVellNat <- function(Radius, RhoParticle, rhoMatrix, DynVisc){
  f_SettlingVelocity (RadCOL, RhoCOL, rhoMatrix, DynVisc,  
                                Matrix = "water")
}

## NOT SURE WHY WE NEED THIS VARIABLE! (Needed for COL and SUSP)