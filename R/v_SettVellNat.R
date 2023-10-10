#'@title SettVellNat
#'@name SettVellNat
#'@description Settling velocity for natural suspended particles
#'@param NaturalRad
#'@param NaturalRho
#'@param rhoMatrix
#'@param DynVisc
#'@export
SettVellNat <- function(RadCOL, RhoCOL, rhoMatrix, DynVisc){
  f_SettlingVelocity (RadCOL, RhoCOL, rhoMatrix, DynVisc,  
                                Matrix = "water")
}
