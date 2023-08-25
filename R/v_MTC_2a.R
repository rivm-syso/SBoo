#'@title PARTIAL MASS TRANSFER COEFFICIENT water and soil to air
#'@name MTC_2a
#'@description Partial mass transfer coefficient from water and soil to air
#'@param WINDspeed Windspeed in compartment/scale [m.s-1]
#'@param MW Molecular weight of compound [g.mol-1]
#'@param from.Matrix Matrix/compartment from which the relevant process is taking place
#'@return 
#'@export
MTC_2a <- function(WINDspeed, MW, Q.10, Temp, T25, kdeg.soil, Matrix){
  ret <- switch(Matrix,
    "water" = 0.01*(0.0004+0.00004*WINDspeed^2)*((0.032/MW)^(0.5*0.5)),
    "soil" = 0.1*Tempfactor.wsds(Q.10,Temp,T25)*kdeg.soil, #KDEG still needs further implementation
    NA
  )
}
