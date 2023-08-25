#'@title Temperature correction degradation rate water/sed/soil
#'@name Tempfactor.wsds
#'@param Q.10 Degradation rate increase factor per /U+00B0 10 C [-]
#'@param Temp Average precipitation [K]
#'@return  Temperature correction degradation rate water/sed/soil [-]
#'@export

Tempfactor.wsds <- function(Q.10,Temp, T25) {

  Q.10^((Temp-T25)/10)
  
}