#' @title General Advection process
#' @name k_Advection
#' @description Calculation of k, given a Flow
#' @param flow advection air
#' @param Volume in m3
#' @return Rate constant for 1rst order process associated with fluxes
#' @export
k_Advection <- function(flow, Volume) {
  flow/Volume
}