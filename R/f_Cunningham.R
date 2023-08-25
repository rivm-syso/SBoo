#' @title Cunningham
#' @name fCunningham
#' @description Collission frequency of ENPs with other particulates due to Gravitational or Intertial impaction in s-1 for heteroagglomeration
#' @param radSpecies radius of the species
#' @return Cunningham 
#' @export
f_Cunningham <- function(radSpecies){
  Knudsen <- (66*10^-9)/(radSpecies)
  1+Knudsen*(1.142+0.558*exp(-0.999/Knudsen))
}
