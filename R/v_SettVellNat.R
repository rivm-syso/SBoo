#'@title SettVellNat
#'@name SettVellNat
#'@description Settling velocity for natural suspended particles
#'@param NaturalRad
#'@param NaturalRho
#'@param rhoMatrix
#'@param DynVisc
#'@export
SettVellNat <- function(RadSPM, RhoSPM, rhoMatrix, DynVisc){
  f_SettlingVelocity (RadSPM, RhoSPM, rhoMatrix, DynVisc,  
                                Matrix = "water")
}
