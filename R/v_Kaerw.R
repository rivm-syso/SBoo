#' @title Aerosol water - air Partition coefficient
#' @name Kaerw
#' @description pm
#' @param Kacompw 
#' @param FRorig
#' @param SubCompartName
#' @return FRears[]
#' @export
Kaerw <- function (Kacompw, FRorig, SubCompartName) {

  switch(SubCompartName,
         "air" = 1/(Kacompw*FRorig),
         NA)
  

  
}