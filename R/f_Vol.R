#' @title Calculate volume of sphere assuming (20181102)
#' @name fVol
#' @description Calculate the volume of spherical particle using the radius in m3
#' @param rad_particle Radius of particle [m]
#' @return fVol [m3]
#' @export
# calculating particle volume:
fVol <- function(rad_particle){ #RADIUS, not diameter!
  (4/3)*pi*rad_particle^3
  # Possible to in future include fractal dimension.
}
