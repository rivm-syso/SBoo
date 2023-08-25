#' @title Rate constant for wet deposition of particulate species due to rain
#' @name k_WetDeposition
#' @description Calculation of the first order rate constant for particle deposition from cloudwater compartment to the soil or water surface [s-1]
#' @param to.AREA Surface area of receiving land or water compartment [m2]
#' @param from.VOLUME Volume of cloudwater compartment [m3]
#'
#'

k_WetDeposition <- function (to.Area, from.Volume, RAINrate)

  (RAINrate* to.Area)/from.Volume

# #Test of function
#to.Area =2.49e-3 * 2.3e11
#from.volume = 3e-7 * 2.3e14
#RAINrate = 2.22e-8

