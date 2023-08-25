#' @title Calculate residence time
#' @name f_TAU
#' @description general approximation of mean residence time
#' @param Area in m2
#' @param WINDspeed in m s-1
#' @return TAU [s]
#' @export
f_TAU <- function (Area, WINDspeed){
  1.5*(0.5 * sqrt(Area*pi/4)/WINDspeed)
}
